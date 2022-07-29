#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/level.h"
#include "../include/cell.h"

/* initialize levels */
enum response init_levels(struct grid *grid)
{
    level * curr_level = grid->first_level;
    // initialize levels of the grid startin from the second level
    for (int id = grid->first_level->id + 1; id <= grid->settings->num_levels; id++)
    {
        if((curr_level->next != NULL))
            exit_with_failure("Error in level.c: Something went wrong, next level should be NULL!");

        // initialize next level
        curr_level->next = (struct level *)malloc(sizeof(struct level));
        if (curr_level->next == NULL)
            exit_with_failure("Error in main.c: Could not allocate memory for next level.");

        level * new_level = curr_level->next;

        new_level->id = id;
        new_level->is_root = false;
        new_level->is_first = false;
        // if its the leaf level (bottom level)
        if(new_level->id == grid->settings->num_levels)
            new_level->is_leaf = true;
        else
            new_level->is_leaf = false;

        new_level->num_cells = (unsigned int)abs(pow(2, grid->settings->num_pivots * id)); // num_cells = 2 ^ (|P| * id)
        new_level->cells = malloc(sizeof(struct cell) * new_level->num_cells);
        if (new_level->cells == NULL)
            exit_with_failure("Error in cell.c: Could not allocate memory for new level cells.");

        for (int i = 0; i < new_level->num_cells; i++)
        {
            new_level->cells[i].id = i;
            new_level->cells[i].center = malloc(sizeof(struct vector));
            new_level->cells[i].center->values = malloc(sizeof(v_type) * grid->settings->num_pivots);
            if (new_level->cells[i].center->values == NULL)
                exit_with_failure("Error in level.c: Could not allocate memory for level center vectors.\n");
        }
        
        // cell edge length = (V / num_cells) ^ 1/|P|
        new_level->cell_edge_length = pow((grid->settings->pivot_space_volume / new_level->num_cells), (1.0/grid->settings->num_pivots));
        // printf("\n¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨");
        // printf("\nNew Level %u, num_cells = %d, cell_edge_length = %f\n", id, new_level->num_cells, new_level->cell_edge_length);    
        
        // (!) link cells in curr_level to cells in new level (curr_cell = current cell in next level)
        for(int i = 0, curr_cell = 0; i < curr_level->num_cells; i++)
        {
            // printf("\n-> Making %d child cells for cell %d in level %d...\n", curr_level->cells[0].num_child_cells, i, curr_level->id);
            // get child cell for current parent cell
            if(&curr_level->cells[i] == NULL)
                exit_with_failure("Error in level.c: Fatal, current cell in current level is a null point.");
            
            struct cell * temp_child_cells = get_child_cells(&curr_level->cells[i], curr_level->cells[i].num_child_cells, new_level, grid->settings);
            if(temp_child_cells == NULL)
                exit_with_failure("Error in level.c: NULL pointer to child cells.");
            // // link child cells to next level
            curr_level->cells[i].children = &new_level->cells[curr_cell];
            for(int j = 0; j < curr_level->cells[i].num_child_cells; j++, curr_cell++)
            {
                // printf("++ cell %d isleaf: %s\n", j, grid->root->cells[0].children[j].is_leaf ? "true" : "false");
                cell_cpy(&new_level->cells[curr_cell], &temp_child_cells[j], grid->settings->num_pivots);
                // print_vector(new_level->cells[curr_cell].center, grid->settings->num_pivots);

                // create filename if its a leaf cell
                if(new_level->is_leaf)
                    create_cell_filename(grid->settings, &new_level->cells[curr_cell]);
            }
            // free memory from temp child cells
            for(int c = curr_level->cells[i].num_child_cells  - 1; c >= 0; c--)
            {
                free(temp_child_cells[c].center->values);
            }
            free(temp_child_cells->center);
            free(temp_child_cells);
        }
        
        // printf("\n¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨");
        new_level->next = NULL;
        curr_level = curr_level->next; // change current last level in list of levels.
    }
    return OK;
}

/* initialize leaf level */
enum response init_first_level(struct grid *grid)
{
    if(grid->root == NULL)
        exit_with_failure("Error in level.c: NULL pointer to root level please run init_root_level() first.");
    
    int i, j;
    /* allocate memory for level */
    grid->first_level = (struct level *)malloc(sizeof(struct level));
    if (grid->first_level == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for level.");

    // link first level and root level
    grid->root->next = grid->first_level;

    grid->first_level->id = 1;
    grid->first_level->num_cells = (unsigned int)abs(pow(2, grid->settings->num_pivots)); // num_cells = 2 ^ (P * 1)
    grid->first_level->cell_edge_length = pow((grid->settings->pivot_space_volume / grid->first_level->num_cells), (1.0/grid->settings->num_pivots));
    grid->first_level->next = NULL;
    
    grid->first_level->is_root = false;
    grid->first_level->is_first = true;
    grid->first_level->is_leaf = false;

    grid->first_level->cells = malloc(sizeof(struct cell) * grid->first_level->num_cells);
    if (grid->first_level->cells == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for first level cells.");

    // printf("\n¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨");
    // printf("\nLevel %u, num_cells = %d, cell_edge_length = %f\n", grid->first_level->id, grid->first_level->num_cells, grid->first_level->cell_edge_length); // cell edge length = (V / num_cells) ^ 1/|P|

    /* create cells in first level */
    for(int c = 0; c < grid->first_level->num_cells; c++)
    {
        init_cell(&grid->first_level->cells[c], grid->first_level->cell_edge_length, grid->first_level->num_cells);
    }

    // allocate memory for center vectors
    vector *center_vectors = malloc(sizeof(struct vector) * grid->first_level->num_cells);
    if (center_vectors == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for list of center vectors.\n");

    for (i = 0; i < grid->first_level->num_cells; i++)
    {
        center_vectors[i].values = malloc(sizeof(v_type) * grid->settings->num_pivots);
        if (center_vectors[i].values == NULL)
            exit_with_failure("Error in level.c: Could not allocate memory for list of center vectors.\n");
        
        // link center vector to cell in first_level
        grid->first_level->cells[i].center = &center_vectors[i];
        grid->first_level->cells[i].level = grid->first_level;
    }

    // link cells in first level to parent cell in root level
    grid->root->cells[0].children = grid->first_level->cells;

    /*
        ndc = number of distinct cooridantes = number of child cells per parent cell = number of cells in the first level 
        compute coordiante for cell center vectors 
        find all distinct coordinate values (v1, v2, v3, ...)
    */
    int ndc = grid->settings->pivot_space_extremity->values[0] / grid->first_level->cell_edge_length; //ndc = number of distinct coordinate values of cell center vectors
    v_type * distinct_coordinates = (v_type *) malloc(sizeof(v_type) * ndc);
    
    if(distinct_coordinates == NULL)
    {
        printf("Error in level.c: Couldn't allocate memory for distinct coordinates, ndc = %d.", ndc);
        exit(-1);
    }
    
    for (i = 0, j = 1; i < ndc; i++, j += 2)
    {
        distinct_coordinates[i] = grid->settings->pivot_space_extremity->values[0] - ((j * grid->first_level->cell_edge_length) / 2);
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
    temp.values = malloc(sizeof(v_type) * grid->settings->num_pivots);
    if (temp.values == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for temp vector.\n");

    int curr_vector = 0;
    create_center_vectors(distinct_coordinates, ndc, grid->settings->num_pivots,
                        grid->settings->num_pivots, center_vectors, temp, 0, &curr_vector);

    // printf("\n-> Making child cells for root cell\n");
    // for(int j = 0; j < grid->root->cells[0].num_child_cells; j++)
    // {
        // printf("++ cell %d isleaf: %s\n", j, grid->root->cells[0].children[j].is_leaf ? "true" : "false");
        // print_vector(grid->root->cells[0].children[j].center, grid->settings->num_pivots);
    // }
    free(temp.values);
    free(distinct_coordinates);
    return OK;
}

/* initialize root level */
enum response init_root(struct grid *grid)
{
    /* allocate memory for level */
    grid->root = (struct level *)malloc(sizeof(struct level));
    if (grid->root == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for root.");

    level * root = grid->root;
    root->id = 0;
    root->num_cells = 1; // root level has one cell
    root->cell_edge_length = pow((grid->settings->pivot_space_volume / root->num_cells), (1.0/grid->settings->num_pivots));
    // root->cell_edge_length = grid->settings->pivot_space_extremity->values[0];
    root->next = NULL;
    
    root->is_root = true;
    root->is_first = false;
    root->is_leaf = false;

    unsigned int num_child_cells = (unsigned int)abs(pow(2, grid->settings->num_pivots));
    root->cells = malloc(sizeof(struct cell) * root->num_cells);
    if (root->cells == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for first level cells.");

    // printf("\n¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨");
    // printf("\nLevel %u, num_cells = %d, cell_edge_length = %f\n", root->id, root->num_cells, root->cell_edge_length); // cell edge length = (V / num_cells) ^ 1/|P|

    /* create cells in root level */
    init_cell(&root->cells[0], root->cell_edge_length, num_child_cells);
    root->cells[0].level = grid->root;
    root->cells[0].center = malloc(sizeof(struct vector));
    root->cells[0].center->values = malloc(sizeof(v_type) * grid->settings->num_pivots);
    
    for(int i = 0; i < grid->settings->num_pivots; i++)
    {
        root->cells[0].center->values[i] = grid->settings->pivot_space_extremity->values[0]/2;
    }

    return OK;
}