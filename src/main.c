#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_loader.h"
#include "../include/gsl_matrix.h"
#include "../include/select_pivots.h"
#include "../include/file_buffer.h"
#include "../include/file_buffer_manager.h"
#include "../include/query_engine.h"

vector * get_dataset_extremity(vector * dataset, unsigned int num_vectors, unsigned int num_dim);
int main()
{
    /* dataset */
    const char *root_directory = "/home/jaouhara/Projects/pexeso-debbug/pexeso/index/";
    const char * bin_files_directory = "/home/jaouhara/Projects/Dissertation/dssdl/encode/binary_files/"; //target directory
    unsigned long num_files = 0ul;
    unsigned int base = 32; // 32 bits to store numbers in binary files
    unsigned int mtr_vector_length = 3, num_dim_metric_space = 3;
    unsigned long long total_vectors = 0ull; // number of vectors in the whole data lake
    unsigned int max_leaf_size = 76; // max vectors in one leaf cell
    double buffered_memory_size = 30; // memory  allocated for file buffers (in MB)
    unsigned int track_vector = 1; // track vectors id (table_id, column_id)

    /* mode 0 = index dataset, 1 = query dataset */
    unsigned int mode = 0;

    /* index settings */
    unsigned int num_levels = 3;  // m
    unsigned int num_pivots = 2;  // number of pivots
    unsigned int fft_scale = 13;   // constant for finding |P| * fft_scale candidate pivots, a good choice of fft_scale is approximately 30 (in paper) and 13 with experiments.

    /* query settings */
    unsigned int join_threshold = 1; // T 
    v_type dist_threshold = 0; // tau


    /* read all vectors in the data set */
    printf("Reading dataset info...");
    get_dataset_info(bin_files_directory, &num_files, &total_vectors, &mtr_vector_length);
    printf("(OK)\n");

    printf("\tNumber of (.bin) files = %lu\n\tNumber of vectors = %llu\n\tVector length in mtric spaces = %u\n\n", num_files, total_vectors, mtr_vector_length);

    printf("Loading dataset files...");
    vector * dataset = load_binary_files(bin_files_directory, 
                                        num_files, total_vectors, base, mtr_vector_length);
        
    if(dataset == NULL)
        exit_with_failure("Error in main.c: Something went wrong, couldn't read dataset vectors!");
    printf("(OK)\n");

    printf("\n\nLooking for pivot vectors...");
    /* search for pivot vectors using pca based algorithm (waiting for response from authors) */
    int dataset_dim [] = {total_vectors, num_dim_metric_space};
    int pivots_mtr_dim [] = {num_pivots, num_dim_metric_space};

    // get pivot vector in from metric dataset
    vector * pivots_mtr = select_pivots(dataset, dataset_dim, num_pivots, fft_scale);
    
    // free dataset
    for(int dv = total_vectors - 1; dv >=0; dv--)
        free(dataset[dv].values);
    free(dataset);  

    printf("(OK)\n");


    /* Transforming data set from metric to pivot space (create distance matrix) */
    printf("\n\nTransforming dataset to pivot space (compute the distance matrix)... ");
    // vector * dataset_ps = map_to_pivot_space(dataset, dataset_dim, pivots_mtr, num_pivots);

    /* map pivot vectors to pivot space pi --> pi' */
    vector * pivots_ps = map_to_pivot_space(pivots_mtr, pivots_mtr_dim, pivots_mtr, num_pivots);
    

    printf("(OK)\n");

    /* pivot space extremity */
    printf("\nextremity vector (in pivot space):\n");
    // vector * pivot_space_extremity = get_rand_vector(num_pivots);
    vector * pivot_space_extremity = get_extremity(pivots_ps, num_pivots);
    print_vector(pivot_space_extremity, num_pivots);
    

    /* initialize index */
    printf("\n\nInitialize index... ");
    struct query_settings * query_settings = init_query_settings(dist_threshold, join_threshold);

    pexeso_index * index = (struct pexeso_index *) malloc(sizeof(struct pexeso_index));
    if (index == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for index!");

    if (!init_index(root_directory, num_pivots, pivots_mtr, pivots_ps, pivot_space_extremity, 
                    num_levels, total_vectors, base, mtr_vector_length, 
                    buffered_memory_size, max_leaf_size, track_vector, 
                    query_settings, index))
        exit_with_failure("Error in main.c: Couldn't initialize index!");

    printf("(OK)\n");

    /* Display settings */
    printf("\n\t\t*** \tINDEX SETTINGS\t ***\t\n");
    printf("--------------------------------------------------------------\n");
    printf("\t\tNumber of pivots = %d\n", index->settings->num_pivots);
    printf("\t\tNumber of levels = %d\n", index->settings->num_levels);
    printf("\t\tPivot space volume = %f\n", index->settings->pivot_space_volume);
    printf("\t\tNumber of leaf cells = %d\n", index->settings->num_leaf_cells);
    printf("\t\tLeaf cell edge length = %f\n", index->settings->leaf_cell_edge_length);
    printf("\t\tPivot vectors (in metric space):\n");
    for(int i = 0; i < num_pivots; i++)
    {
        print_vector(&index->settings->pivots_mtr[i], num_dim_metric_space);
    }
    printf("\t\tPivot vectors (in pivot space):\n");
    for(int i = 0; i < num_pivots; i++)
    {
        print_vector(&index->settings->pivots_ps[i], num_pivots);
    }
    printf("--------------------------------------------------------------\n");


    /* Build levels */
    printf("\n\nBuild levels... ");
    if (!init_first_level(index))
        exit_with_failure("Error in main.c: Couldn't initialize first level!");

    if(index->first_level == NULL)
        exit_with_failure("Error in main.c: first level not initialized for index!");    

    if (!init_levels(index))
        exit_with_failure("Error in main.c: Couldn't initialize index levels!");
    printf("(OK)\n");

    
    /* insert dataset in index */
    /* read all vectors in the data set */
    printf("Index dataset vectors...");
    if (!index_binary_files(index, bin_files_directory, num_files, base))
        exit_with_failure("Error in main.c: Couldn't initialize first level!");
    printf("(OK)\n");

    /* print index */
    print_index(index);

    /* write index to disk */
    if (!index_write(index))
        exit_with_failure("Error main.c:  Could not save the index to disk.\n");
    

    /* destroy index */
    if (!index_destroy(index, index->first_level))
        exit_with_failure("Error main.c:  Could not destroy index.\n");
    
    // exit(1);

    // free index and index settings
    for(int p = index->settings->num_pivots - 1; p >= 0; p--)
        free(index->settings->pivots_mtr[p].values);
    free(index->settings->pivots_mtr);

    for(int p = index->settings->num_pivots - 1; p >= 0; p--)
        free(index->settings->pivots_ps[p].values);
    free(index->settings->pivots_ps);

    free(index->settings->pivot_space_extremity->values);
    free(index->settings->pivot_space_extremity);

    free(index->settings->query_settings);
    free(index->settings);
    free(index);

    return 0;
}

vector * get_dataset_extremity(vector * dataset, unsigned int num_vectors, unsigned int num_dim)
{
    vector * extremity = malloc(sizeof(struct vector));
    extremity->values = malloc(sizeof(v_type)*num_dim);
    for(int j = 0; j < num_dim; j++)
    {
        extremity->values[j] = FLT_MIN;
    }

    for(int v = 0; v < num_vectors; v++)
    {
        for(int j = 0; j < num_dim; j++)
        {
            if(dataset[v].values[j] > extremity->values[j])
                extremity->values[j] = dataset[v].values[j];
        }
        
    }

    return extremity;

}