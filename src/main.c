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


    /* mode 0 = index dataset, 1 = query dataset */
    unsigned int mode = 0;

    /* index settings */
    unsigned int num_levels = 2; // m
    unsigned int num_pivots = 2; // number of pivots

    /* get set of pivot vectors*/
    vector * pivots = (struct vector *) malloc(sizeof(struct vector) * num_pivots);
    for (int i = 0; i < num_pivots; i++)
    {
        pivots[i].values = (float *)malloc(sizeof(v_type) * num_pivots);
    }

    // search for pivot vectors using pca based algorithm. 
    // (...)

    vector * pivot_space_extrimity = get_extrimity(pivots, num_pivots);


    /* initialize index */
    pexeso_index * index = (struct pexeso_index *) malloc(sizeof(struct pexeso_index));
    if (index == NULL)
        exit_with_error("Error in main.c: Couldn't allocate memory for index!");

    if (!init_index(root_directory, num_pivots, pivot_space_extrimity, num_levels, index))
        exit_with_error("Error in main.c: Couldn't initialize index!");
    


    

    printf("Number of leaf cells = %d\n", index->settings->num_leaf_cells);
    printf("Leaf cell edge length = %.2f\n", index->settings->leaf_cell_edge_length);
    printf("Index has been initialized!\n");

    if (!init_first_level(index))
        exit_with_error("Error in main.c: Couldn't initialize first level!");

    if(index->first_level == NULL)
        exit_with_error("Error in main.c: first level not initialized for index!");
    
    if (!init_levels(index))
        exit_with_error("Error in main.c: Couldn't initialize index levels!");
    

    /* get image of the data lake in the pivot space */
    // metric_to_pivot_space(char * dataset_dir)







    // if (!index_binary_files(index, bin_files_directory, num_files, base))
    //     exit_with_error("Error in main.c: Couldn't initialize first level!");

    return 0;
}

/* print error and kill process */
void exit_with_error(char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}
