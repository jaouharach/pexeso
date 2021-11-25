#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"

/* initialize levels */
response init_levels(pexeso_index *index)
{
    level *level = index->first_level;
    // initialize levels of the index startin from the second level
    for (int id = 2; id <= index->settings->num_levels; id++)
    {
        if((level->next != NULL))
            exit_with_error("Error in level.c: Something went wrong, next level should be NULL!");

        // initialize next level
        level->next = (struct level *)malloc(sizeof(struct level));

        if (level->next == NULL)
            exit_with_error("Error in main.c: Could not allocate memory for next level.");

        level->next->id = id;
        // num_cells = 2 ^ (P * id)
        level->next->num_cells = (unsigned int)abs(pow(2, index->settings->num_dim * id));

        printf("Num cells in level %u = %d\n", id, level->next->num_cells);

        level->next->cell_edge_length = (index->settings->max_coordinate - index->settings->min_coordinate) / pow(level->next->num_cells, (1 / index->settings->num_dim));
        level->next->next = NULL;
        level->next->prev = level; // point back to last level in index

        level = level->next; // change current last level in list of levels.
    }

    // initialize cells not leaf cells.
    // init_leaf_cells(level->next, index->settings);

    return OK;
}

/* initialize leaf level */
response init_first_level(pexeso_index *index)
{
    // initialize first level
    index->first_level = (struct level *)malloc(sizeof(struct level));

    if (index->first_level == NULL)
        exit_with_error("Error in main.c: Could not allocate memory for level.");

    index->first_level->id = 1;
    // num_cells = 2 ^ (P * 1)
    index->first_level->num_cells = (unsigned int)abs(pow(2, index->settings->num_dim));

    printf("Num cells in first level = %d\n", index->first_level->num_cells);

    index->first_level->cell_edge_length = (index->settings->max_coordinate - index->settings->min_coordinate) / pow(index->first_level->num_cells, (1 / index->settings->num_dim));
    index->first_level->next = NULL;
    index->first_level->prev = NULL;

    // initialize cells.
    // init_leaf_cells(index->first_level, index->settings);

    return OK;
}
