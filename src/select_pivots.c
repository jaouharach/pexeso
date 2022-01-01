#include "../include/index.h"
#include "../include/select_pivots.h"

/* pivot selection algorithm */
vector * select_pivots(vector * dataset, int * dataset_dim, unsigned int num_pivots, unsigned int fft_scale)
{
    // run fft to get a candidate set of outliers
    int num_cp = num_pivots*fft_scale; // number of candidate pivots

    // get a set of candidate pivots (outliers)
    vector * candidate_pivots = fft(dataset, dataset_dim, num_cp);

    // compute the distance matrix (num_vector x num_cp) to the current set of candidate pivots (outliers)
    vector * dataset_ps = map_to_pivot_space(dataset, dataset_dim, candidate_pivots, num_cp);
    
    return candidate_pivots;
}


/* FFT: Farthest First Traversal, k = number of outliers to be found */
vector *fft(vector *data_set, int * dataset_dim, unsigned int k)
{
    int num_vectors = dataset_dim[0]; // total vector in the dataset
    int v_len = dataset_dim[1]; // vector length

    vector *outliers = malloc(sizeof(struct vector) * k);
    if(outliers == NULL)
        exit_with_failure("Error in select_pivots.c: Couldn't allocate memory for outliers.");

    for (int i = 0; i < k; i++)
    {
        outliers[i].values = calloc(v_len, sizeof(v_type));
        if(outliers[i].values == NULL)
            exit_with_failure("Error in select_pivots.c: Couldn't allocate memory for outlier values.");

    }

    float bsf_d;
    int bsf_v;

    if (k >= num_vectors)
        warning("Warning in globals.c: number of outliers is greater or equal to total vectors in the dataset.");

    // pick an arbitrary vector as first outlier
    unsigned int rand_idx = (rand() % (num_vectors - 1));
    vector_cpy(&outliers[0], &data_set[rand_idx], v_len);

    // find k outliers, skip first
    for (int i = 1; i < k; i++)
    {
        bsf_d = FLT_MIN; // max distance to the current set of outliers
        bsf_v = -1;      // position of the farthest vector in the data set.
        // find the farthest vector to the vectors in outliers list
        for (int j = 0; j < num_vectors; j++)
        {
            // check if vector is already in oultiers
            // (...)
            float d = min_distance(outliers, i + 1, &data_set[j], v_len);
            if (d > bsf_d)
            {
                bsf_d = d;
                bsf_v = j;
            }
        }
        outliers[i].values = data_set[bsf_v].values;
    }
    return outliers;
}
// min distance to a set of outliers (for fft)
float min_distance(vector *outliers, unsigned int num_outliers, vector *v, unsigned int v_len)
{
    float min_d = FLT_MAX, d;
    for (int i = 0; i < num_outliers; i++)
    {
        d = euclidean_distance(v, &outliers[i], v_len);
        if (d < min_d)
            min_d = d;
    }
    return min_d;
}

/* map vector to pivot space  vector --> vector_mapping*/
void map_vector(vector *v, unsigned int v_len, vector *v_mapping, vector *pivots, unsigned int num_pivots)
{
    for (int i = 0; i < num_pivots; i++)
        v_mapping->values[i] = euclidean_distance(v, &pivots[i], v_len);   
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
        
        map_vector(&dataset_mtr[i], num_dim_metric_space, &dataset_ps[i], pivots, num_pivots);
    }
    return dataset_ps;
}
