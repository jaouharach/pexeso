#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_loader.h"
#include "../include/select_pivots.h"

vector * get_dataset_extremity(vector * dataset, unsigned int num_vectors, unsigned int num_dim);
int main()
{
    /* dataset */
    const char *root_directory = "/home/jaouhara/Projects/pexeso/index";
    const char * bin_files_directory = "/home/jaouhara/Projects/Dissertation/dssdl/encode/binary_files/"; //target directory
    unsigned int num_files = 5;
    unsigned int base = 32; // 32 bits to store numbers in binary files
    unsigned int num_dim_metric_space = 3;
    unsigned long total_vectors = 0ul; // number of vectors in the whole data lake

    /* mode 0 = index dataset, 1 = query dataset */
    unsigned int mode = 0;

    /* index settings */
    unsigned int num_levels = 2;  // m
    unsigned int num_pivots = 2;  // number of pivots
    unsigned int fft_scale = 1;   // constant for finding |P| * c candidate pivots

    /* read all vectors in the data set */
    printf("Reading data set vectors...");
    vector * dataset = load_binary_files(&total_vectors, bin_files_directory, 
                                        num_files, base, num_dim_metric_space);
        
    if(dataset == NULL)
        exit_with_failure("Error in main.c: Couldn't read dataset vectors!.");
    printf("(OK)\n");
    printf("Total vectors = %lu\n", total_vectors);


    printf("\n\nLooking for pivot vectors... ");
    /* search for pivot vectors using pca based algorithm (waiting for response from authors) */
    // temp solution: use fft (fathest fist traversal to find k  pivots)
    int dataset_dim [] = {total_vectors, num_dim_metric_space};
    vector * pivots_mtr = select_pivots(dataset, dataset_dim, num_pivots, fft_scale);
    printf("(OK)\n");

    /* Mapping pivot vectors from metric to pivot space */
    printf("\n\nMapping pivots vectors... ");
    vector * pivots = malloc(sizeof(struct vector) * num_pivots);
    for(int i = 0; i < num_pivots; i++)
    {
        pivots[i].values = malloc(sizeof(v_type) * num_pivots);
        if(pivots[i].values == NULL)
            exit_with_failure("Error in main.c: couldn't allocate memory for pivot mapping.");
        transform_vector(&pivots_mtr[i], num_dim_metric_space,
                        &pivots[i], pivots_mtr, num_pivots);   
    }
    printf("(OK)\n");

    
    printf("Pivots found using FFT:\n");
    printf("--> In metric space:\n");
    for(int j = 0; j < num_pivots; j++)
        print_vector(&pivots_mtr[j], num_dim_metric_space);
    printf("\n\n");
    printf("--> In pivot space:\n");
    for(int j = 0; j < num_pivots; j++)
        print_vector(&pivots[j], num_pivots);

    /* pivot space extremity */
    printf("extremity vector in pivot space:\n");
    vector * pivot_space_extremity = get_extremity(pivots, num_pivots);
    print_vector(pivot_space_extremity, num_pivots);
    

    /* Transforming data set from metric to pivot space (create distance matrix) */
    printf("\n\nTransforming dataset to pivot space (compute the distance matrix)... ");
    vector * dataset_img = malloc(sizeof(struct vector) * total_vectors);
    for(int i = 0; i < total_vectors; i++)
    {
        dataset_img[i].values = malloc(sizeof(v_type) * num_pivots);
        transform_vector(&dataset[i], num_dim_metric_space,
                        &dataset_img[i], pivots_mtr, num_pivots);
        // printf("V%d:\n", i+1);
        // print_vector(&dataset[i], num_dim_metric_space);
        // printf("V%d*:\n", i+1);
        // print_vector(&dataset_img[i], num_pivots);
    }
    printf("(OK)\n");

    /* find extremity in dataset */

    // printf("\n\nextremity vector in pivot space:\n");
    // pivot_space_extremity = get_dataset_extremity(dataset_img, total_vectors, num_pivots);
    // print_vector(pivot_space_extremity, num_pivots);
    // exit(1);

    /* initialize index */
    printf("\n\nInitialize index... ");
    pexeso_index * index = (struct pexeso_index *) malloc(sizeof(struct pexeso_index));
    if (index == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for index!");

    if (!init_index(root_directory, num_pivots, pivot_space_extremity, num_levels, index))
        exit_with_failure("Error in main.c: Couldn't initialize index!");
    printf("(OK)\n");

    /* Display settings */
    printf("\n\t\t*** \tINDEX SETTINGS\t ***\t\n");
    printf("--------------------------------------------------------------\n");
    printf("\t\tNumber of pivots = %d\n", index->settings->num_dim);
    printf("\t\tNumber of levels = %d\n", index->settings->num_levels);
    printf("\t\tPivot space volume = %f\n", index->settings->pivot_space_volume);
    printf("\t\tNumber of leaf cells = %d\n", index->settings->num_leaf_cells);
    printf("\t\tLeaf cell edge length = %f\n", index->settings->leaf_cell_edge_length);
    printf("--------------------------------------------------------------\n");


    /* BUild levels */
    printf("\n\nBuild levels... ");
    if (!init_first_level(index))
        exit_with_failure("Error in main.c: Couldn't initialize first level!");

    if(index->first_level == NULL)
        exit_with_failure("Error in main.c: first level not initialized for index!");    

    /* (!) implement code to create other levels and  link parent cells in first level to child cells in other levels. */
    
    if (!init_levels(index))
        exit_with_failure("Error in main.c: Couldn't initialize index levels!");

    printf("(OK)\n");

    exit(1);
    // if (!index_binary_files(index, bin_files_directory, num_files, base))
    //     exit_with_failure("Error in main.c: Couldn't initialize first level!");

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