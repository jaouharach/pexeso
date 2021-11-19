#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"


/* initialize leaf level */
int init_leaf_level(index_settings *settings, level *level)
{
    //num_cells = (max - min) ^ num_dim /cell_length ^ num_dim
    
    // check if num cells is integer
    if (fmod(
            pow(settings->max_coordinate - settings->min_coordinate, settings->num_dim),
            pow(settings->leaf_cell_edge_length, settings->num_dim) == 0))
    {
        level->id = 0;
        level->num_cells = (unsigned int)abs(
            pow(settings->max_coordinate - settings->min_coordinate, settings->num_dim) / pow(settings->leaf_cell_edge_length, settings->num_dim));
        
        printf("num cells: %d\n", level->num_cells);

        level->cell_length = settings->leaf_cell_edge_length;
        level->next_level = NULL;

        // initialize cells.
        init_leaf_cells(level, settings);

        return OK;
    }
    else
        exit_with_error("Warning: Number of cells can only be integer! please change settings.");
}

