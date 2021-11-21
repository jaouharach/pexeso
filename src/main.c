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
    const char *root_directory = "/home/jaouhara/Projects/pexeso/index";
    unsigned int num_dim = 3;
    float max_coordinate = 6; // we assume than all dims have the same min and max coordinate
    float min_coordinate = 0;
    float leaf_cell_edge_length = 2;
    int curr_state = OK;
    unsigned int mode = 0;
    const char * bin_files_directory = "/home/jaouhara/Projects/Dissertation/dssdl/encode/binary_files/";
    unsigned int num_files = 5;
    unsigned int base = 32; // 32 bits to store numbers in binary files
    
    // initialize index
    pexeso_index * index = (struct pexeso_index *) malloc(sizeof(struct pexeso_index));

    if (index == NULL)
        exit_with_error("Error in main.c: Couldn't allocate memory for index!");

    if (!init_index(root_directory, num_dim, max_coordinate, min_coordinate, leaf_cell_edge_length, index))
        exit_with_error("Error in main.c: Couldn't initialize index!");
  
    printf("Index has been initialized!\n");

    if (!init_leaf_level(index))
        exit_with_error("Error in main.c: Couldn't initialize first level!");

    if(index->first_level == NULL)
        exit_with_error("Error in main.c: first level not initialized for index!");
        
    if (!index_binary_files(index, bin_files_directory, num_files, base))
        exit_with_error("Error in main.c: Couldn't initialize first level!");

    return 0;
}

/* print error and kill process */
void exit_with_error(char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}
