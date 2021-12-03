#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"

/* initialize levels */
response init_levels(pexeso_index *index)
{
    printf("Init leves...\n");
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
        // num_cells = 2 ^ (|P| * id)
        level->next->num_cells = (unsigned int)abs(pow(2, index->settings->num_dim * id));
        level->next->cells = (struct cell *)malloc(sizeof(struct cell) * level->next->num_cells);
        if (level->next->cells == NULL)
            exit_with_error("Error in main.c: Could not allocate memory for next level cells.");

        // cell edge length = (V / num_cells) ^ 1/|P|
        level->next->cell_edge_length = pow((index->settings->pivot_space_volume / level->next->num_cells), (1/index->settings->num_dim));
        printf("Num cells in level %u = %d, cell edge length = %.2f\n", id, level->next->num_cells, level->next->cell_edge_length);    

        level->next->next = NULL;
        level->next->prev = level; // point back to last level in index

        // initialize cells.
        for(int c = 0; c < level->next->num_cells; c++)
        {
            init_cell(&level->next->cells[c], level->next->cell_edge_length);
        }
        level = level->next; // change current last level in list of levels.
    }

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
    index->first_level->cells = (struct cell *)malloc(sizeof(struct cell) * index->first_level->num_cells);
        if (index->first_level->cells == NULL)
            exit_with_error("Error in main.c: Could not allocate memory for first level cells.");

    // cell edge length = (V / num_cells) ^ 1/|P|
    index->first_level->cell_edge_length = pow((index->settings->pivot_space_volume / index->first_level->num_cells), (1/index->settings->num_dim));
    printf("Num cells in first level = %d, cell edge length = %.2f\n", index->first_level->num_cells, index->first_level->cell_edge_length);    

    index->first_level->next = NULL;
    index->first_level->prev = NULL;

    // initialize cells.
    for(int c = 0; c < index->first_level->num_cells; c++)
    {
        init_cell(&index->first_level->cells[c], index->first_level->cell_edge_length);
    }
    return OK;
}
