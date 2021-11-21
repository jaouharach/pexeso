#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"


/* initialize leaf level */
response init_leaf_level(pexeso_index * index)
{
    // initialize first level
    level * level = (struct level *)malloc(sizeof(struct level));

    if (level == NULL)
        exit_with_error("Error in main.c: Could not allocate memory for level.");

    //num_cells = (max - min) ^ num_dim /cell_edge_length ^ num_dim
    
    // check if num cells is integer
    if (fmod(
            pow(index->settings->max_coordinate - index->settings->min_coordinate, index->settings->num_dim),
            pow(index->settings->leaf_cell_edge_length, index->settings->num_dim) == 0))
    {
        level->id = 0;
        level->num_cells = (unsigned int)abs(
            pow(index->settings->max_coordinate - index->settings->min_coordinate, index->settings->num_dim) / pow(index->settings->leaf_cell_edge_length, index->settings->num_dim));
        
        printf("num cells: %d\n", level->num_cells);

        level->cell_edge_length = index->settings->leaf_cell_edge_length;
        level->next = NULL;
        level->prev = NULL;
        
        // initialize cells.
        init_leaf_cells(level, index->settings);

        index->first_level = level;

        return OK;
    }
    else
        exit_with_error("Error in level.c: Number of cells can only be integer! please change settings.");
}

