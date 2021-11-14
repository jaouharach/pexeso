#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "headers/main.h"

int main()
{
    const char *root_directory = "/home/jaouhara/Projects/pexeso/index";
    unsigned int num_dim = 2;
    float max_dim_value = 1; // we assume than all dims have the same min and max value
    float min_dim_value = -1;
    float leaf_cell_length = 1;
    int curr_state = OK;

    // initialize index
    index *pexeso_index = (index *)malloc(sizeof(struct index));

    if (pexeso_index == NULL)
        exit_with_error("Couldn't allocate memory for index!");

    if (!init_index(root_directory, num_dim, max_dim_value, min_dim_value, leaf_cell_length, pexeso_index))
        exit_with_error("Warning: Couldn't initialize index!");

    // printf("Root directory = %s\n", pexeso_index->settings->root_directory);
    // printf("Vector size = %d\n", pexeso_index->settings->num_dim);

    printf("Index has been initialized!\n");

    // initialize first level
    level *level = (struct level *)malloc(sizeof(struct level));
    if (!init_leaf_level(pexeso_index->settings, level))
        exit_with_error("Warning: Couldn't initialize first level!");

    printf("Leaf level has been initialized!\n");

    return 0;
}

int init_leaf_cells(level * leaf_level, unsigned int num_dim)
{
    leaf_level->cells = (cell *)malloc(sizeof(struct cell) * leaf_level->num_cells);
    int num_row = sqrt(leaf_level->num_cells);

    for (int i = 0; i < num_row; i++)
    {
        for (int j = 0; j < num_row; j++)
        {
            //initialize cell
            init_leaf_cell(&leaf_level->cells[i + j], leaf_level->cell_length);
            leaf_level->cells[i + j].center = (vector *) malloc(sizeof(struct vector));
            leaf_level->cells[i + j].center->values = (float *) malloc(sizeof(float) * num_dim);
                
            //set center vector
            for (int k = 0; k < num_dim; k++)
            {
                leaf_level->cells[i + j].center->values[k] = 0;
            }
        }
    }

    return OK;
}

int init_leaf_cell(cell *cell, float length)
{
    cell->filepath = "";
    cell->parent = NULL;
    cell->is_leaf = 1;
    cell->length = length;
    cell->num_vectors = 0;
    cell->center = NULL;

    return OK;
}

void exit_with_error(char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

int init_leaf_level(index_settings *settings, level *level)
{
    // num_cells = (max - min) ^ num_dim /cell_length ^ num_dim
    // check if num cells is integer
    if (fmod(
            pow(settings->max_dim_value - settings->min_dim_value, settings->num_dim),
            pow(settings->leaf_cell_length, settings->num_dim) == 0))
    {
        level->id = 0;
        level->num_cells = (unsigned int)abs(
            pow(settings->max_dim_value - settings->min_dim_value, settings->num_dim) / pow(settings->leaf_cell_length, settings->num_dim));
        printf("num cells: %d\n", level->num_cells);
        level->cell_length = settings->leaf_cell_length;
        level->next_level = NULL;

        // initialize cells.
        init_leaf_cells(level, settings->num_dim);

        return OK;
    }
    else
        exit_with_error("Warning: Number of cells can only be integer! please change settings.");
}

int init_index(const char *root_directory, unsigned int num_dim,
               float max_dim_value, float min_dim_value,
               float leaf_cell_length, index *pexeso_index)
{
    pexeso_index->settings = (index_settings *)malloc(sizeof(index_settings));
    if (pexeso_index->settings == NULL)
        return FAILED;

    pexeso_index->first_level = 0;
    pexeso_index->total_records = 0;

    pexeso_index->settings->root_directory = root_directory;
    pexeso_index->settings->num_dim = num_dim;
    pexeso_index->settings->max_dim_value = max_dim_value;
    pexeso_index->settings->min_dim_value = min_dim_value;
    pexeso_index->settings->leaf_cell_length = leaf_cell_length;

    return OK;
}