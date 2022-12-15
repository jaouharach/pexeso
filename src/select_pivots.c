#include "../include/hgrid.h"
#include "../include/gsl_matrix.h"
#include "../include/select_pivots.h"
#include "../include/stats.h"
#include <float.h>
#include <time.h>

/* alloc dataset for pivot selection */
// vector * alloc_dataset_memory(struct grid_settings * settings, unsigned long long total_vectors)
// {
//     unsigned long num_bytes = settings->mtr_buffered_memory_size * 1024 * 1024;
//     unsigned long vector_size_in_bytes = sizeof(v_type) * settings->mtr_vector_length;
//     long max_buffered_vectors = (long)(num_bytes / vector_size_in_bytes);

//     // allocate a big chunk of memory
//     char * memory_array = calloc(max_buffered_vectors, vector_size_in_bytes);
//     if (memory_array == NULL)
//         exit_with_failure("Error in selcet_pivots.c:"
//                           " Cannot allocate allocate dataset memory for pivot selection.\n");

//     // split memory into vectors
//     char * current_record = memory_array;
//     unsigned int current_record_index  = 0;
//     vector * dataset = memory_array;
//     unsigned int next = 0;
//     if(current_record_index > max_buffered_vectors)
//         exit_with_failure("Error in select_pivots: Excceded max memory for storing dataset vectors!");

//     for(int v = 0; v < total_vectors; v++)
//     {
//         current_record_index++;
//         for(int i = 0; i < settings->mtr_vector_length; i++)
//         {
            
//         }
//         next = next + (settings->mtr_vector_length *  sizeof(v_type));
//     }
// }
/* select pivots using best fft scale */
struct vector *select_pivots_with_best_fft_scale(struct vector *dataset, int *dataset_dim, int *dims, unsigned int min_fft, unsigned int max_fft, unsigned short num_iter)
{
    unsigned int num_pivots = dims[0];
    unsigned int metric_dim = dims[1];

    struct best_fft bsf;
    bsf.exremity = FLT_MIN;
    bsf.pivots_mtr = malloc(sizeof(struct vector) * num_pivots);
    for(int p = 0; p < num_pivots; p++)
        bsf.pivots_mtr[p].values = malloc(sizeof(v_type) * metric_dim);

    vector * pivots_mtr = NULL;
    vector * pivots_ps = NULL;
    vector * ps_extremity = NULL;

    for(unsigned int f = min_fft; f <= max_fft; f++)
    {
        // printf("FFT SCALE = %u", f); 
        for(int i = 0; i < num_iter; i++)
        {
            pivots_mtr = select_pivots(dataset, dataset_dim, num_pivots, f);
            pivots_ps = map_to_pivot_space(pivots_mtr, dims, pivots_mtr, num_pivots);
            ps_extremity = get_extremity(pivots_ps, num_pivots);

            if (ps_extremity->values[0] > bsf.exremity)
            {
                bsf.exremity = ps_extremity->values[0];
                bsf.fft_scale = f;

                // copy pivots to bsf
                for(int p = 0; p < num_pivots; p++)
                    for(int m = 0; m < metric_dim; m++)
                    {
                        bsf.pivots_mtr[p].values[m] = pivots_mtr[p].values[m];
                    }
            }
            for(int p = 0; p < num_pivots; p++)
            {
                free(pivots_mtr[p].values);
                free(pivots_ps[p].values);
            }
            free(pivots_mtr);
            free(pivots_ps);
            pivots_mtr = NULL;
            pivots_ps = NULL;

            printf("\rFFT = %d, extremity = %.2f (%d/%d)", f, ps_extremity->values[0],  i+1, num_iter);
            fflush(stdout);
            // print_vector(ps_extremity, num_pivots);
            free(ps_extremity->values);
            free(ps_extremity);
        }
    }
    if(pivots_ps != NULL)
    {
        for(int p = 0; p < num_pivots; p++)
        {
            free(pivots_ps[p].values);
        }
        free(pivots_ps);
    }

    printf("\nBest scale = {fft = %u, extremity = %.2f}\n", bsf.fft_scale, bsf.exremity);
    return bsf.pivots_mtr;
}
/* pivot selection algorithm */
vector * select_pivots(vector * dataset, int * dataset_dim, unsigned int num_pivots, unsigned int fft_scale)
{
    // run fft to get a candidate set of outliers
    int num_cp = num_pivots*fft_scale; // number of candidate pivots

    // get a set of candidate pivots (outliers)
    vector * candidate_pivots = fft(dataset, dataset_dim, num_cp);

    // compute the distance matrix (num_vector x num_cp) to the current set of candidate pivots (outliers)
    vector * dataset_ps = map_to_pivot_space(dataset, dataset_dim, candidate_pivots, num_cp);
    
    // compute PCA with EM method
    int dim_pcset [] = {num_pivots, num_cp+1};
    gsl_matrix * pcset = empca(dataset_ps, dataset_dim[0], num_cp, num_pivots);

    // select pivots from the pca result, result = indecies of best pivots in candidate_pivots
    int * result = select_pivots_by_pca_result_angle(pcset, dim_pcset, num_pivots);

    vector * pivots = malloc(sizeof(struct vector) * num_pivots);
    for(int i = 0; i < num_pivots; i++)
    {
        pivots[i].values =  malloc(sizeof(v_type) * dataset_dim[1]);
        // printf("index of P%d = %d\n", i, result[i]);
        for(int j = 0; j < dataset_dim[1]; j++)
        {
            pivots[i].values[j] = candidate_pivots[result[i]].values[j];
        }
    }

    // free memory
    for (int i = num_cp - 1; i >= 0 ; i--)
        free(candidate_pivots[i].values);
    free(candidate_pivots);

    for (int i = dataset_dim[0] - 1 ; i >= 0 ; i--)
        free(dataset_ps[i].values);
    free(dataset_ps);

    
    free(result);
    gsl_matrix_free(pcset);

    return pivots;
}

/* FFT: Farthest First Traversal, k = number of outliers to be found */
vector *fft(vector *data_set, int * dataset_dim, unsigned int k)
{
    int num_vectors = dataset_dim[0]; // total vector in the dataset
    int v_len = dataset_dim[1]; // vector length

    // allocate memory for outliers
    int * outliers_idx = malloc(sizeof(int) * k); // indexes of outiers in the dataset
    struct vector *outliers = malloc(sizeof(struct vector) * k);
    if(outliers == NULL || outliers_idx == NULL)
        exit_with_failure("Error in select_pivots.c: Couldn't allocate memory for outliers.");

    for (int i = 0; i < k; i++)
    {       
        outliers_idx[i] = -1; // to make that this outlier is not found yet
        outliers[i].values = malloc(sizeof(v_type) * v_len);
        if(outliers[i].values == NULL)
            exit_with_failure("Error in select_pivots.c: Couldn't allocate memory for outlier values.");
    }

    float bsf_d;
    int bsf_v;

    if (k >= num_vectors)
        warning("Warning in globals.c: number of outliers is greater or equal to total vectors in the dataset.");

    // pick an arbitrary vector as first outlier
    srand (time(NULL));
    unsigned int rand_idx = (rand() % (num_vectors - 1));
    outliers_idx[0] = rand_idx;
    // printf("\nrand idx = %d\n", rand_idx);
    vector_cpy(&outliers[0], &data_set[rand_idx], v_len);
    
    // printf("first outlier: \n");
    // print_vector(&outliers[0], v_len);


    // find the other  (k-1) outliers, skip first
    for (int i = 1; i < k; i++)
    {
        bsf_d = -FLT_MAX; // max distance to the current set of outliers
        bsf_v = -1;      // position of the farthest vector in the data set.
        // find the farthest vector to the vectors in outliers list
        for (int j = 0; j < num_vectors; j++)
        {   
            // check if vector is already in oultiers
            if(in_outliers(outliers_idx, j, k) == 1)
                continue;

            float d = min_distance(outliers, i, &data_set[j], v_len);

            if(d == FLT_MAX)
            {
                // print_vector(&data_set[j], v_len);
                exit_with_failure("Error int select_pivots.c: min distance cannot be equal to FLT_MAX!");
            }

            if (d > bsf_d) // bsf_d is the fartest distance to the current set of outliers
            {
                bsf_d = d;
                bsf_v = j;
            }
        }
        if(bsf_v == -1)
            exit_with_failure("Error in select_pivots.c: Something went wrong! cannot add outlier of index -1, no such vector in the dataset.");
        
        for(int v = 0; v < v_len; v++)
            outliers[i].values[v] = data_set[bsf_v].values[v];
    }

    // free memory
    free(outliers_idx);
    return outliers;
}

// check if a new outlier is already in the list of outliers
unsigned int in_outliers(int * outliers_idx, int new_outlier, int num_outliers)
{
    for(int o = 0; o < num_outliers; o++)
    {
        if(outliers_idx[o] == new_outlier)
            return 1;
    }
    return 0;
}

// min distance to a set of outliers (for fft)
float min_distance(vector *outliers, unsigned int num_outliers, vector *v, unsigned int v_len)
{
    float min_d = FLT_MAX, d;
    for (int i = 0; i < num_outliers; i++)
    {
        d = euclidean_distance_cmp(v, &outliers[i], v_len);
        if (d < min_d)
        {
            min_d = d;
        }    
    }
    return min_d;
}

/* map vector to pivot space  v --> v'*/
enum response map_vector(struct vector *v, unsigned int v_len, struct vector *v_mapping, struct grid_settings * settings, bool in_query_grid)
{
    float d;
    int out_of_pivot_space = 0;
    struct vector * pivots_mtr = settings->pivots_mtr;
    struct vector * extremity = settings->pivot_space_extremity;
    unsigned int num_pivots = settings->num_pivots;
    for (int i = 0; i < num_pivots; i++)
    {
        d = euclidean_distance(v, &pivots_mtr[i], v_len);
        if(d > extremity->values[i])
            out_of_pivot_space = 1;
        v_mapping->values[i] = d; 
    } 

    if(out_of_pivot_space == 1)
    {
        if(in_query_grid)
            COUNT_NEW_OUT_OF_PIVOT_SPACE_QUERY_VECTOR
        else
            COUNT_NEW_OUT_OF_PIVOT_SPACE_VECTOR
    }

    return OK;
}

/*  get distance matrix of a dataset (mapping of vectors to the pivot space), 
    dataset_mtr (dataset in metric space) --> dataset_ps (dataset in pivot space) */
vector * map_to_pivot_space(vector * dataset_mtr, int * dataset_dim, vector * pivots, unsigned int num_pivots)
{
    int num_vectors = dataset_dim[0];
    int num_dim_metric_space = dataset_dim[1];

    // allocate memory for dataset in pivot space
    vector * dataset_ps = malloc(sizeof(struct vector) * num_vectors);
    if(dataset_ps == NULL)
        exit_with_failure("Error in select_pivots.c: Couldn't allocate memory for new dataset in pivot space.");
    for(int i = 0; i < num_vectors; i++)
    {
        dataset_ps[i].values = malloc(sizeof(v_type) * num_pivots);
        if(dataset_ps[i].values == NULL)
            exit_with_failure("Error in select_pivots.c: couldn't allocate memory for values of vector mapping.");
        
        // map to pivot space
        for (int j = 0; j < num_pivots; j++)
        {
            dataset_ps[i].values[j] = euclidean_distance(&dataset_mtr[i], &pivots[j], num_dim_metric_space);
        } 
    }
    return dataset_ps;
}

/* EMPCA: Expectation Management for Principal Component Analysis, return pc matrix */
gsl_matrix * empca(vector *data_set, unsigned int num_vectors, unsigned int dim, unsigned int num_pc)
{
    unsigned int num_iteration = 20;
    
    int dim_data[] = {num_vectors, dim};
    int dim_data_t[] = {dim, num_vectors};

    // centerize data
    gsl_matrix * data = to_gsl_matrix(data_set, dim_data);
    gsl_matrix_centerize(data, dim_data);

    // transpose matrix of the data_set
    gsl_matrix * data_t = gsl_matrix_get_transpose(data, dim_data);

    // printf("Mean centerised data matrix:\n");
    // gsl_matrix_print(data, dim_data);

    // printf("transpose of data matrix:\n");
    // gsl_matrix_print(data_t, dim_data_t);

    // create a (dim x num_pc) matrix of random weigths, each weigth wij = [-0.5, 0.5]
    int dim_c[] = {dim, num_pc};
    int dim_ct[] = {num_pc, dim};
    int dim_ct_mult_c[] = {num_pc, num_pc};

    gsl_matrix * C = gsl_matrix_alloc(dim, num_pc);
    gsl_matrix * Ct = gsl_matrix_alloc(num_pc, dim);

    if(C == NULL || Ct == NULL)
        exit_with_failure("Error in select_pivots.c: Couldn't allocate memory for empca matrix.");

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < num_pc; j++)
            gsl_matrix_set(C, i, j, ((float)rand() / (float)((unsigned)RAND_MAX + 1)) - 0.5);  
    

    int dim_x [] = {num_pc, num_vectors};
    int dim_xt [] = {num_vectors, num_pc};

    gsl_matrix *x = NULL;
    gsl_matrix *xt = gsl_matrix_alloc(num_vectors, num_pc);

    // printf("C matrix:\n");
    // gsl_matrix_print(C, dim_c);

    // repeat 20 times
    for (int i = 0; i < num_iteration; i++)
    {
        // compute Ct
        if(Ct != NULL) gsl_matrix_free(Ct);
        Ct = gsl_matrix_get_transpose(C, dim_c);        

        // compute x, x = [(Ct * C) ^-1 * Ct] * data
        if(x != NULL) gsl_matrix_free(x);

        // (Ct * C) ^-1
        gsl_matrix * Ct_C_inv = gsl_matrix_inverse(gsl_matrix_product(Ct, dim_ct, C, dim_c), dim_ct_mult_c);
        // [(Ct * C) ^-1 * Ct]
        gsl_matrix * Ct_C_inv_Ct = gsl_matrix_product(Ct_C_inv, dim_ct_mult_c, Ct, dim_ct);
        // [(Ct * C) ^-1 * Ct] * data
        x = gsl_matrix_product(Ct_C_inv_Ct,dim_ct, data_t, dim_data_t);

        // free memory
        gsl_matrix_free(Ct_C_inv);
        gsl_matrix_free(Ct_C_inv_Ct);

        // compute xt 
        if(xt != NULL) gsl_matrix_free(xt);

        xt = gsl_matrix_get_transpose(x, dim_x);

        // compute C, C = (data * Xt) * (X * Xt) ^-1
        if(C != NULL) gsl_matrix_free(C);
        
        // (X * Xt) ^-1
        gsl_matrix * X_Xt_inv = gsl_matrix_inverse(gsl_matrix_product(x, dim_x, xt, dim_xt), dim_ct_mult_c);
        // (data * Xt)
        gsl_matrix * data_Xt = gsl_matrix_product(data_t, dim_data_t, xt, dim_xt);
         // (data * Xt) * (X * Xt) ^-1
        C = gsl_matrix_product(data_Xt, dim_c, X_Xt_inv, dim_ct_mult_c);

        // free memory
        gsl_matrix_free(X_Xt_inv);
        gsl_matrix_free(data_Xt);
    }
    

    // orthonormalization of matrix C (SVD)
    orthonormalization(C, dim_c, 1);

    // cov_matrix = Covariance matrix of (Ct * data)t
    int dim_ctd[] = {num_pc, num_vectors}; // dim of (Ct * data)
    int dim_ctdt [] = {num_vectors, num_pc}; // dim of (Ct * data)t
    // int dim_cov [] = {num_pc, num_pc};
    int dim_v [] = {num_pc, num_pc};

    // compute Ct
    if(Ct != NULL) gsl_matrix_free(Ct);
    Ct = gsl_matrix_get_transpose(C, dim_c);

    // compute (Ct * data)t   
    gsl_matrix * ct_mult_data =  gsl_matrix_product(Ct, dim_ct, data_t, dim_data_t);
    gsl_matrix * ct_mult_data_t = gsl_matrix_get_transpose(ct_mult_data, dim_ctd);

    // compute covariance matrix
    gsl_matrix * cov_matrix = gsl_matrix_covariance(ct_mult_data_t, dim_ctdt);

    // Eigen value decomposition of cov_matrix: cov_matrix = V * [eigen_vals] * Vt
    gsl_vector *eigen_vals = gsl_vector_alloc(num_pc);
    gsl_matrix *V = gsl_matrix_alloc(num_pc, num_pc);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(num_pc);

    gsl_eigen_symmv(cov_matrix, eigen_vals, V, w);
    

    // result
    int dim_result [] = {dim+1, num_pc};
    gsl_matrix * result = gsl_matrix_alloc(dim+1, num_pc);

    // (C x Vf)
    gsl_matrix * Vf = gsl_matrix_flip_columns(V, dim_v);
    gsl_matrix * c_mult_vf = gsl_matrix_alloc(dim, num_pc);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C, Vf, 0.0, c_mult_vf); 

    // assign (C * Vf) to sub matrix (dim x num_pc) starting at (1,0)
    gsl_matrix_set_part(result, 1, 0, c_mult_vf, dim, num_pc);

    // assign eigen values for first row in result matrix
    gsl_vector * eigen_vals_fliped = vector_flip(eigen_vals, num_pc);
    gsl_matrix_set_row(result, 0,  eigen_vals_fliped);

    gsl_matrix * result_t = gsl_matrix_get_transpose(result, dim_result);

    // free memory
    gsl_matrix_free(data);
    gsl_matrix_free(data_t);
    gsl_matrix_free(C);
    gsl_matrix_free(Ct);
    gsl_matrix_free(ct_mult_data);
    gsl_matrix_free(x);
    gsl_matrix_free(xt);
    gsl_matrix_free(V);
    gsl_matrix_free(Vf);
    gsl_matrix_free(cov_matrix);
    gsl_matrix_free(ct_mult_data_t);
    gsl_vector_free(eigen_vals_fliped);
    gsl_vector_free(eigen_vals);
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(result);

    return result_t;
}

/* select pivots from the empca result (pc with highest projection on dataset)
    returns indecies of the pivots with the highest projection */
int * select_pivots_by_pca_result_angle(gsl_matrix * pca_result, int *pca_result_dim, int num_pivots)
{
    unsigned int num_pc = pca_result_dim[0];
    unsigned int num_cols = pca_result_dim[1] - 1;

    int dim_pc [] = {num_pc+1, num_cols};
    int dim_pc_t [] = {num_cols, num_pc+1};
    gsl_matrix * PC = gsl_matrix_alloc(num_pc+1, num_cols);

    gsl_matrix_set_part(PC, 1, 0, gsl_matrix_get_part(pca_result, 0, 1, num_pc, num_cols), num_pc, num_cols);

    // convert matrix to positive matrix
    gsl_to_positive_matrix(PC, dim_pc);
    // extra column to sort changes on pct matrix
    for (int i=0; i < num_cols; i++)
            gsl_matrix_set(PC, 0, i, (double)i);
    
    // printf("PCt matrix:\n");
    // gsl_matrix_print(gsl_matrix_get_transpose(PC, dim_pc), dim_pc_t);
    
    int counter = 0;
    int * result = calloc(num_pivots, sizeof(*result));
    if(result == NULL)
        exit_with_failure("Error in select_pivots.c: Couldn't allocate memory for integer array.");
    
    gsl_matrix * pc_view = NULL;
    gsl_matrix * PCt = gsl_matrix_get_transpose(PC, dim_pc);

    //for jth pc, find the ith axis d(., pi), on which the pc has the largest projection
    for(int i = 0; (i < num_cols) && (counter < num_pivots); i++)
    {
        for(int j = 1; (j <= num_pc) && (counter < num_pivots); j++)
            {
                // PCt
                pc_view = gsl_matrix_sort_by_column(PCt, dim_pc_t, j);
                // printf("PCt matrix sorted by column %d:\n", j);
                // gsl_matrix_print(pc_view, dim_pc_t);
                
                int point = (int) gsl_matrix_get(pc_view, num_cols -1, 0);
                
                if(array_add(result, counter, point))
                {   
                    // printf("highest projection in row = %d\n", point);
                    counter ++;
                }
                // // replace max value in column j with zero
                // gsl_matrix_set(pc_view, num_cols -1, j, 0);
                gsl_matrix_free(pc_view);
            }
    }

    // free memory
    gsl_matrix_free(PC);
    gsl_matrix_free(PCt);

    return result; // indecies of axis d(., pi) of which pcs have the highest projection 
}
