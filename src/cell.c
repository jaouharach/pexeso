#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"

/* initialize leaf cells in a level */
response init_leaf_cells(level *leaf_level, index_settings *settings)
{
    /* allocate memory for cells */
    leaf_level->cells = (cell *)malloc(sizeof(struct cell) * leaf_level->num_cells);
    int i, j, ndc = log(leaf_level->num_cells) / log(settings->num_dim); //ndc = number of distinct coordinate values of cell center vectors

    /* copmute distinct coordinate values (v1, v2, v3, ...) */
    float *distinct_coordinates = (float *)malloc(sizeof(float) * ndc);
    for (i = 0, j = 1; i < ndc; i++, j += 2)
    {
        distinct_coordinates[i] = settings->max_coordinate - ((j * settings->leaf_cell_edge_length) / 2);
    }

    /* init cells with their center vectors */
    vector temp;
    temp.values = malloc(sizeof(float) * settings->num_dim);

    if (temp.values == NULL)
        exit_with_error("Error in cell.c: Could not allocate memory for temp vector.\n");

    vector *center_vectors = malloc(sizeof(struct vector) * leaf_level->num_cells);

    if (center_vectors == NULL)
        exit_with_error("Error in cell.c: Could not allocate memory for list of center vectors.\n");
    

    for (i = 0; i < leaf_level->num_cells; i++)
    {
        center_vectors[i].values = malloc(sizeof(float) * settings->num_dim);

        if (center_vectors[i].values == NULL)
            exit_with_error("Error in cell.c: Could not allocate memory for list of center vectors.\n");
        
    }

    init_center_vectors(distinct_coordinates, ndc, settings->num_dim,
                        settings->num_dim, center_vectors, temp, 0);

    for (i = 0; i < leaf_level->num_cells; i++)
    {
        init_leaf_cell(&leaf_level->cells[i], leaf_level->cell_length);
        leaf_level->cells[i].center = &center_vectors[i];

        /* printf("(");
        for(j = 0; j < settings->num_dim; j++)
            printf("%.2f, ", leaf_level->cells[i].center->values[j]);
        printf(")\n"); */
    }

    return OK;
}

/* initialize leaf cell */
response init_leaf_cell(cell *cell, float length)
{
    cell->filename = "";
    cell->parent = NULL;
    cell->is_leaf = 1;
    cell->edge_length = length;
    cell->num_vectors = 0;
    cell->center = NULL;

    return OK;
}

/* create center vectors for all cells of the data space */
void init_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector *center_vectors, vector temp, int append_at)
{
    static int curr_vector = 0;
    if (ndc == 0 || k > ndc)
    {
        exit_with_error("Error in cell.c: Couldn'l initialize center vector!\n");
    }

    if (k == 0)
    {
        // add temp to sub_array
        for (int x = 0; x < dim; x++)
            center_vectors[curr_vector].values[x] = temp.values[x];
        curr_vector++;
        return;
    }

    for (int j = 0; j < ndc; j++)
    {
        // add distinct_coordinates[j] to temp values
        temp.values[append_at] = distinct_coordinates[j];
        init_center_vectors(distinct_coordinates, ndc, k - 1, dim, center_vectors, temp, append_at + 1);
    }
}