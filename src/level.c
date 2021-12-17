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
            exit_with_failure("Error in level.c: Something went wrong, next level should be NULL!");

        // initialize next level
        level->next = (struct level *)malloc(sizeof(struct level));

        if (level->next == NULL)
            exit_with_failure("Error in main.c: Could not allocate memory for next level.");

        level->next->id = id;
        // num_cells = 2 ^ (|P| * id)
        level->next->num_cells = (unsigned int)abs(pow(2, index->settings->num_dim * id));
        
        // cell edge length = (V / num_cells) ^ 1/|P|
        level->next->cell_edge_length = pow((index->settings->pivot_space_volume / level->next->num_cells), (1.0/index->settings->num_dim));
        printf("\n**************************************************************");
        printf("\nLevel %u, num_cells = %d, cell_edge_length = %f\n", id, level->next->num_cells, level->next->cell_edge_length);    

        make_level_cells(level->next, index->settings);

        level->next->next = NULL;
        level->next->prev = level; // point back to last level in index
        level = level->next; // change current last level in list of levels.
    }
    return OK;
}

/* initialize leaf level */
response init_first_level(pexeso_index *index)
{
    /* allocate memory for level */
    index->first_level = (struct level *)malloc(sizeof(struct level));
    if (index->first_level == NULL)
        exit_with_failure("Error in main.c: Could not allocate memory for level.");

    index->first_level->cells = (struct cell *)malloc(sizeof(struct cell) * index->first_level->num_cells);
    if (index->first_level->cells == NULL)
        exit_with_failure("Error in main.c: Could not allocate memory for first level cells.");


    index->first_level->id = 1;
    index->first_level->num_cells = (unsigned int)abs(pow(2, index->settings->num_dim)); // num_cells = 2 ^ (P * 1)
    index->first_level->cell_edge_length = pow((index->settings->pivot_space_volume / index->first_level->num_cells), (1.0/index->settings->num_dim));
    index->first_level->next = NULL;
    index->first_level->prev = NULL;
    
    printf("\n**************************************************************");
    printf("\nLevel 0, num_cells = %d, cell_edge_length = %f\n", index->first_level->num_cells, index->first_level->cell_edge_length); // cell edge length = (V / num_cells) ^ 1/|P|

    /* initialize cells*/
    for(int c = 0; c < index->first_level->num_cells; c++)
    {
        init_cell(&index->first_level->cells[c], index->first_level->cell_edge_length);
    }

    /* compute coordiante for cell center vectors */
    int i, j, ndc = index->settings->pivot_space_extrimity->values[0] / index->first_level->cell_edge_length; //ndc = number of distinct coordinate values of cell center vectors

    // find all distinct coordinate values (v1, v2, v3, ...)
    v_type * distinct_coordinates = (v_type *) malloc(sizeof(v_type) * ndc);
    for (i = 0, j = 1; i < ndc; i++, j += 2)
    {
        distinct_coordinates[i] = index->settings->pivot_space_extrimity->values[0] - ((j * index->first_level->cell_edge_length) / 2);
    }
    printf("number of distinct coordinates = %d\n", ndc);
    printf("Distinct coordinates:\n");
    printf("(");
    for(i = 0; i < ndc; i++)
    {
        printf("%f, ", distinct_coordinates[i]);
    }
    printf(")\n");

    /* create center vectors */
    vector temp;
    temp.values = malloc(sizeof(v_type) * index->settings->num_dim);
    if (temp.values == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for temp vector.\n");

    // allocate memory for center vectors
    vector *center_vectors = malloc(sizeof(struct vector) * index->first_level->num_cells);
    if (center_vectors == NULL)
        exit_with_failure("Error in level.c: Could not allocate memory for list of center vectors.\n");

    for (i = 0; i < index->first_level->num_cells; i++)
    {
        center_vectors[i].values = malloc(sizeof(v_type) * index->settings->num_dim);
        if (center_vectors[i].values == NULL)
            exit_with_failure("Error in level.c: Could not allocate memory for list of center vectors.\n");
        
        index->first_level->cells[i].center = &center_vectors[i];
    }

    init_center_vectors(distinct_coordinates, ndc, index->settings->num_dim,
                        index->settings->num_dim, center_vectors, temp, 0);

    // print center vectors
    // for (i = 0; i < index->first_level->num_cells; i++)
    // {
    //     printf("(");
    //     for(j = 0; j < index->settings->num_dim; j++)
    //         printf("%.2f, ", index->first_level->cells[i].center->values[j]);
    //     printf(")\n");
    // }

    return OK;
}
