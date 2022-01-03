#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"

/* initialize levels */
response init_levels(pexeso_index *index)
{
    level * curr_level = index->first_level;
    // initialize levels of the index startin from the second level
    for (int id = 2; id <= index->settings->num_levels; id++)
    {
        if((curr_level->next != NULL))
            exit_with_failure("Error in level.c: Something went wrong, next level should be NULL!");

        // initialize next level
        curr_level->next = (struct level *)malloc(sizeof(struct level));
        if (curr_level->next == NULL)
            exit_with_failure("Error in main.c: Could not allocate memory for next level.");

        level * new_level = curr_level->next;

        new_level->id = id;
        new_level->num_cells = (unsigned int)abs(pow(2, index->settings->num_pivots * id)); // num_cells = 2 ^ (|P| * id)
        new_level->cells = malloc(sizeof(struct cell) * new_level->num_cells);
        if (new_level->cells == NULL)
            exit_with_failure("Error in cell.c: Could not allocate memory for new level cells.");

        // allocate memory for center vectors
        vector *center_vectors = malloc(sizeof(struct vector) * new_level->num_cells);
        if (center_vectors == NULL)
            exit_with_failure("Error in level.c: Could not allocate memory for list of center vectors.\n");

        for (int i = 0; i < new_level->num_cells; i++)
        {
            center_vectors[i].values = malloc(sizeof(v_type) * index->settings->num_pivots);
            if (center_vectors[i].values == NULL)
                exit_with_failure("Error in level.c: Could not allocate memory for list of center vectors.\n");
            
            new_level->cells[i].center = &center_vectors[i];
        }
        // cell edge length = (V / num_cells) ^ 1/|P|
        new_level->cell_edge_length = pow((index->settings->pivot_space_volume / new_level->num_cells), (1.0/index->settings->num_pivots));
        printf("\n¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨");
        // printf("\nCurr Level %u, num_cells = %d, cell_edge_length = %f\n", curr_level->id, curr_level->num_cells, curr_level->cell_edge_length);    
        printf("\nNew Level %u, num_cells = %d, cell_edge_length = %f\n", id, new_level->num_cells, new_level->cell_edge_length);    
        
        // (!) link cells in curr_level to cells in new level (curr_cell = current cell in next level)
        for(int i = 0, curr_cell = 0; i < curr_level->num_cells; i++)
        {
            printf("\n-> Making %d child cells for cell %d in level %d...\n", curr_level->cells[0].num_child_cells, i, curr_level->id);
            // get child cell for current parent cell
            if(&curr_level->cells[i] == NULL)
                exit_with_failure("Error in level.c: Fatal, current cell in current level is a null point.");
            cell * child_cells = get_child_cells(&curr_level->cells[i], curr_level->cells[i].num_child_cells, index->settings);
            
            // // link child cells to next level
            curr_level->cells[i].children = &new_level->cells[curr_cell];
            for(int j = 0; j < curr_level->cells[i].num_child_cells; j++, curr_cell++)
            {
                printf("++ Adding cell %d, center vector: \n", curr_cell);
                cell_cpy(&new_level->cells[curr_cell], &child_cells[j], index->settings->num_pivots);
                print_vector(new_level->cells[curr_cell].center, index->settings->num_pivots);
            }
        }
        printf("\n¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨");
        new_level->next = NULL;
        new_level->prev = curr_level; // point back to last level in index
        curr_level = curr_level->next; // change current last level in list of levels.
    }
    return OK;
}

/* initialize leaf level */
response init_first_level(pexeso_index *index)
{
    int i, j;

    /* allocate memory for level */
    index->first_level = (struct level *)malloc(sizeof(struct level));
    if (index->first_level == NULL)
        exit_with_failure("Error in main.c: Could not allocate memory for level.");

    index->first_level->id = 1;
    index->first_level->num_cells = (unsigned int)abs(pow(2, index->settings->num_pivots)); // num_cells = 2 ^ (P * 1)
    index->first_level->cell_edge_length = pow((index->settings->pivot_space_volume / index->first_level->num_cells), (1.0/index->settings->num_pivots));
    index->first_level->next = NULL;
    index->first_level->prev = NULL;
    
    index->first_level->cells = malloc(sizeof(struct cell) * index->first_level->num_cells);
    if (index->first_level->cells == NULL)
        exit_with_failure("Error in main.c: Could not allocate memory for first level cells.");

    printf("\n¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨");
    printf("\nLevel 1, num_cells = %d, cell_edge_length = %f\n", index->first_level->num_cells, index->first_level->cell_edge_length); // cell edge length = (V / num_cells) ^ 1/|P|

    /* create cells in first level */
    for(int c = 0; c < index->first_level->num_cells; c++)
    {
        init_cell(&index->first_level->cells[c], index->first_level->cell_edge_length, index->first_level->num_cells);
    }

    // allocate memory for center vectors
    vector *center_vectors = malloc(sizeof(struct vector) * index->first_level->num_cells);
    if (center_vectors == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for list of center vectors.\n");

    for (i = 0; i < index->first_level->num_cells; i++)
    {
        center_vectors[i].values = malloc(sizeof(v_type) * index->settings->num_pivots);
        if (center_vectors[i].values == NULL)
            exit_with_failure("Error in level.c: Could not allocate memory for list of center vectors.\n");
        
        // link center vector to cell in first_level
        index->first_level->cells[i].center = &center_vectors[i];
    }

    /*
        ndc = number of distinct cooridantes = number of child cells per parent cell = number of cells in the first level 
        compute coordiante for cell center vectors 
        find all distinct coordinate values (v1, v2, v3, ...)
    */
    int ndc = index->settings->pivot_space_extremity->values[0] / index->first_level->cell_edge_length; //ndc = number of distinct coordinate values of cell center vectors
    v_type * distinct_coordinates = (v_type *) malloc(sizeof(v_type) * ndc);
    if(distinct_coordinates == NULL)
        exit_with_failure("Error in level.c: Couldn't allocate memory for distinct coordinates.");

    for (i = 0, j = 1; i < ndc; i++, j += 2)
    {
        distinct_coordinates[i] = index->settings->pivot_space_extremity->values[0] - ((j * index->first_level->cell_edge_length) / 2);
    }
    
    // printf("number of distinct coordinates = %d\n", ndc);
    // printf("Distinct coordinates:\n");
    // printf("(");
    // for(i = 0; i < ndc; i++)
    // {
    //     printf("%f, ", distinct_coordinates[i]);
    // }
    // printf(")\n");

    /* create center vectors */
    vector temp;
    temp.values = malloc(sizeof(v_type) * index->settings->num_pivots);
    if (temp.values == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for temp vector.\n");

    create_center_vectors(distinct_coordinates, ndc, index->settings->num_pivots,
                        index->settings->num_pivots, center_vectors, temp, 0);

    return OK;
}
