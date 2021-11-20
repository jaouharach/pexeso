#include <stdio.h>
#include <stdlib.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_buffer.h"

/* initialize index */
response init_index(const char *root_directory, unsigned int num_dim,
                    float max_coordinate, float min_coordinate,
                    float leaf_cell_edge_length, index *index)
{
    index->settings = (index_settings *)malloc(sizeof(index_settings));
    if (index->settings == NULL)
        exit_with_error("Error in index.c: Couldn't allocate memory for index settings!");

    index->first_level = 0;
    index->total_records = 0;

    index->settings->root_directory = root_directory;
    index->settings->num_dim = num_dim;
    index->settings->max_coordinate = max_coordinate;
    index->settings->min_coordinate = min_coordinate;
    index->settings->leaf_cell_edge_length = leaf_cell_edge_length;

    return OK;
}

/* append vector to index */
response append_vector(index * index, vector * vector)
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
        cell =
            cell_route_to_closest_child(cell, vector, index->settings->num_dim);

        if(cell == NULL)
            exit_with_error("Error in index.c: Could not route to closest child cell.\n");
    }

    // add vector to leaf cell file buffer.
    int s = cell->file_buffer->buffered_list_size;
    if (s == 0)
    {
        cell->file_buffer->buffered_list = NULL;
        cell->file_buffer->buffered_list = malloc(sizeof(struct vector *));

        if (cell->file_buffer->buffered_list == NULL)
            exit_with_error("Error in index.c: Could not"
                            "allocate memory for the buffered list.\n");
    }
    else
    {
        cell->file_buffer->buffered_list = realloc(cell->file_buffer->buffered_list,
                                                   sizeof(struct vector *) * cell->file_buffer->buffered_list_size + 1);
        if (cell->file_buffer->buffered_list == NULL)
            exit_with_error("Error in index.c: Could not"
                            "reallocate memory for the buffered list.\n");
    }
    
    cell->file_buffer->buffered_list[s].table_id = vector->table_id;
    cell->file_buffer->buffered_list[s].set_id = vector->set_id;

    for (int i = 0; i < index->settings->num_dim; ++i)
    {
        cell->file_buffer->buffered_list[s].values[i] = vector->values[i];
    }

    cell->file_buffer->buffered_list_size++;

    return OK;
}

response index_binary_files(char * dataset_dir, unsigned int l, index * index)
{

    return OK;
}
