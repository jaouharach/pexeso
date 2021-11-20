#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"

int main()
{
    const char *root_directory = "/home/jaouhara/Projects/pexeso/index";
    unsigned int num_dim = 3;
    float max_coordinate = 6; // we assume than all dims have the same min and max coordinate
    float min_coordinate = 0;
    float leaf_cell_edge_length = 2;
    int curr_state = OK;
    unsigned int mode = 0;

    // initialize index
    index * index = (struct index *)malloc(sizeof(struct index));

    if (index == NULL)
        exit_with_error("Error in main.c: Couldn't allocate memory for index!");

    if (!init_index(root_directory, num_dim, max_coordinate, min_coordinate, leaf_cell_edge_length, index))
        exit_with_error("Error in main.c: Couldn't initialize index!");
  
    printf("Index has been initialized!\n");

    // initialize first level
    level * level = (struct level *)malloc(sizeof(struct level));

    if (level == NULL)
        exit_with_error("Error in main.c: Could not allocate memory for level.\n");

    if (!init_leaf_level(index->settings, level))
        exit_with_error("Error in main.c: Couldn't initialize first level!");

    printf("Leaf level has been initialized!\n");

    return 0;
}

/* print error and kill process */
void exit_with_error(char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}
