#include <stdio.h>
#include <stdlib.h>
#include "../include/index.h"

/* initialize index */
int init_index(const char *root_directory, unsigned int num_dim,
               float max_coordinate, float min_coordinate,
               float leaf_cell_edge_length, index *pexeso_index)
{
    pexeso_index->settings = (index_settings *)malloc(sizeof(index_settings));
    if (pexeso_index->settings == NULL)
        return FAILED;

    pexeso_index->first_level = 0;
    pexeso_index->total_records = 0;

    pexeso_index->settings->root_directory = root_directory;
    pexeso_index->settings->num_dim = num_dim;
    pexeso_index->settings->max_coordinate = max_coordinate;
    pexeso_index->settings->min_coordinate = min_coordinate;
    pexeso_index->settings->leaf_cell_edge_length = leaf_cell_edge_length;

    return OK;
}

