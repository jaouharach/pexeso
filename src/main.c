#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_loader.h"

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
    unsigned int num_levels = 2; // m
    unsigned int num_pivots = 2; // number of pivots


    /* read all vectors in the data set */
    printf("Reading data set vectors...\n");
    vector * data_set = load_binary_files(&total_vectors, bin_files_directory, 
                                        num_files, base, num_dim_metric_space);
        
    if(data_set == NULL)
        exit_with_failure("Error in main.c: Couldn't read dataset vectors!.");

    printf("Finish reading data set vectors, total vectors = %lu\n", total_vectors);

    /* get set of pivot vectors*/
    // search for pivot vectors using pca based algorithm (waiting for response from authors)
    /* 
        ...............................
        ...............................
        ...............................
        ...............................
        ...............................
        ...............................
        ...............................
        ...............................
    */

    // temp solution: use fft (fathest fist traversal to find k  pivots)
    vector * pivots = fft(data_set, total_vectors, num_pivots, num_dim_metric_space);
    printf("\n\nPivots found using FFT:\n");
    for(int j = 0; j < num_pivots; j++)
        print_vector(&pivots[j], num_dim_metric_space);

    /* pivot space extrimity */
    printf("\n\nExtrimity vector in pivot space:\n");
    vector * pivot_space_extrimity = get_extrimity(pivots, num_pivots);
    print_vector(pivot_space_extrimity, num_pivots);
    exit(1);

    /* initialize index */
    pexeso_index * index = (struct pexeso_index *) malloc(sizeof(struct pexeso_index));
    if (index == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for index!");

    if (!init_index(root_directory, num_pivots, pivot_space_extrimity, num_levels, index))
        exit_with_failure("Error in main.c: Couldn't initialize index!");
    

    printf("Number of leaf cells = %d\n", index->settings->num_leaf_cells);
    printf("Leaf cell edge length = %.2f\n", index->settings->leaf_cell_edge_length);
    printf("Index has been initialized!\n");

    if (!init_first_level(index))
        exit_with_failure("Error in main.c: Couldn't initialize first level!");

    if(index->first_level == NULL)
        exit_with_failure("Error in main.c: first level not initialized for index!");
    
    if (!init_levels(index))
        exit_with_failure("Error in main.c: Couldn't initialize index levels!");
    

    /* get image of the data lake in the pivot space */
    // metric_to_pivot_space(char * dataset_dir)


    // if (!index_binary_files(index, bin_files_directory, num_files, base))
    //     exit_with_failure("Error in main.c: Couldn't initialize first level!");

    return 0;
}

