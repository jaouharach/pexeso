#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_buffer.h"

/* initialize index */
response init_index(const char * root_directory,
                unsigned int num_pivots,
                vector * extremity,
                unsigned int num_levels,
                pexeso_index * index)
{
    index->settings = (index_settings *)malloc(sizeof(index_settings));
    if (index->settings == NULL)
        exit_with_failure("Error in index.c: Couldn't allocate memory for index settings!");

    index->first_level = NULL;
    index->total_records = 0;

    index->settings->root_directory = root_directory;
    index->settings->num_dim = num_pivots; // number of dimensions is equal to the number of pivots
    index->settings->pivot_space_extremity = extremity;
    index->settings->num_levels = num_levels;
    index->settings->num_leaf_cells = pow(2, num_pivots * num_levels); // 2^(|P| * m) number of cells depends on num_pivots to ensure same length in all edges.

    
    // volume of the pivot space = multiplication of all extremity coordinates.
    index->settings->pivot_space_volume = 1.0;
    for(int i = 0; i < num_pivots; i++)
    {
        index->settings->pivot_space_volume *= fabs(extremity->values[i]);
    }

    // leaf cell edge length = (V / num_leaf_cells) ^ 1/|P|
    index->settings->leaf_cell_edge_length = pow((index->settings->pivot_space_volume / index->settings->num_leaf_cells), (1.0/num_pivots));
    
    return OK;
}
/* get extremity vector of the pivot space, pivots must be outliers for extremity to be accurate */
vector * get_extremity(vector * pivot_vectors, unsigned int num_pivots)
{
    // extremity =  farthest vector in the pivot space (holds maximum distance to a pivot)
    vector * pivot_space_extremity = malloc(sizeof(struct vector));
    pivot_space_extremity->values = (v_type *) malloc(sizeof(v_type) * num_pivots);

    v_type max = FLT_MIN;
    for(int i = 0; i < num_pivots; i++)
    {
        for(int j = 0; j < num_pivots; j++)
            if(pivot_vectors[i].values[j] > max)
                max = pivot_vectors[i].values[j];
    }
    for(int i = 0; i < num_pivots; i++)
        pivot_space_extremity->values[i] = max;
        
    return pivot_space_extremity;
}

/* append vector to index */
response insert_vector(pexeso_index * index, vector *vector)
{
    // find the closest cell in first level.
    float bsf = FLT_MAX;
    cell *cell = NULL;
    for (int c = 0; c < index->first_level->num_cells; c++)
    {
        float d = euclidean_distance(vector, index->first_level->cells[c].center, index->settings->num_dim);
        if (d <= bsf)
        {
            bsf = d;
            cell = &index->first_level->cells[c];
        }
    }

    // loop children of cell untill you find closest leaf cell.
    while (!cell->is_leaf)
    {
        cell = cell_route_to_closest_child(cell, vector, index->settings->num_dim);

        if (cell == NULL)
            exit_with_failure("Error in index.c: Could not route to closest child cell.\n");
    }

    // add vector to leaf cell file buffer.
    // allocate memory for new vector
    int s = cell->file_buffer->buffered_list_size;
    if (s == 0)
    {
        cell->file_buffer->buffered_list = NULL;
        cell->file_buffer->buffered_list = malloc(sizeof(struct vector *));

        if (cell->file_buffer->buffered_list == NULL)
            exit_with_failure("Error in index.c: Could not"
                            "allocate memory for the buffered list.\n");

        cell->file_buffer->buffered_list[s].values = (v_type *) malloc(sizeof(v_type) * index->settings->num_dim);
        
        if (cell->file_buffer->buffered_list[s].values == NULL)
            exit_with_failure("Error in index.c: Could not"
                            "reallocate memory for buffered list.\n");
    }
    else
    {
        /* Resize memory allocate for buffered list to fit for the new vector */
        cell->file_buffer->buffered_list = realloc(cell->file_buffer->buffered_list,
                                                   sizeof(struct vector) * (cell->file_buffer->buffered_list_size + 1));
        if (cell->file_buffer->buffered_list == NULL)
            exit_with_failure("Error in index.c: Could not"
                            "reallocate memory for buffered list.\n");

        cell->file_buffer->buffered_list[s].values = (v_type *) malloc(sizeof(v_type) * index->settings->num_dim);
        
        if (cell->file_buffer->buffered_list[s].values == NULL)
            exit_with_failure("Error in index.c: Could not"
                            "reallocate memory for buffered list.\n");
    }

    // add vector to buffered list
    cell->file_buffer->buffered_list[s].table_id = vector->table_id;
    cell->file_buffer->buffered_list[s].set_id = vector->set_id;
    for (int i = 0; i < index->settings->num_dim; ++i)
    {
        cell->file_buffer->buffered_list[s].values[i] = vector->values[i];
    }

    cell->file_buffer->buffered_list_size++;

    return OK;
}

void display_indexed_vectors(pexeso_index * index)
{
    for(int i = 0; i < index->settings->num_levels; i++)
    {
        
    }
}