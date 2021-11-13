#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "main.h"

int main()
{
    const char *root_directory = "/home/jaouhara/Projects/pexeso/index";
    unsigned int vector_size = 2;
    float max_dim_value = 1; // we assume than all dims have the same min and max value
    float min_dim_value = -1;
    float leaf_cell_area = .5;
    int curr_state = OK;

    // initialize index
    index *pexeso_index = (index *)malloc(sizeof(struct index));

    if (pexeso_index == NULL)
        exit_with_error("Couldn't allocate memory for index!");

    if (!init_index(root_directory, vector_size, max_dim_value, min_dim_value, leaf_cell_area, pexeso_index))
        exit_with_error("Warning: Couldn't initialize index!");

    printf("Index has been initialized!\n");

    // initialize first level
    level *level = (struct level *)malloc(sizeof(struct level));
    if (!init_level(pexeso_index->settings, level))
        exit_with_error("Warning: Couldn't initialize first level!");

    printf("First level has been initialized!\n");

    // printf("Root directory = %s\n", pexeso_index->settings->root_directory);
    // printf("First level = %d\n", pexeso_index->first_level);
    // printf("Vector size = %d\n", pexeso_index->settings->vector_size);

    return 0;
}

int init_level(index_settings *settings, level *level)
{
    // num_cells = area /cell_area ^2
    // check if num cells is integer
    if (fmod(pow(settings->max_dim_value - settings->min_dim_value, 2), settings->leaf_cell_area) == 0)
    {
        level->id = 0;
        level->num_cells = (unsigned int) abs(pow(settings->max_dim_value - settings->min_dim_value, 2) / settings->leaf_cell_area);
        printf("num cells: %d\n", level->num_cells);
        level->cell_area = settings->leaf_cell_area;
        level->cells = NULL;
        level->next_level = NULL;

        return OK;
    }
    else
        exit_with_error("Warning: Number of cells can only be integer! please change settings.");
}

int init_index(const char *root_directory, unsigned int vector_size,
               float max_dim_value, float min_dim_value,
               float leaf_cell_area, index *pexeso_index)
{
    pexeso_index->settings = (index_settings *)malloc(sizeof(index_settings));
    if (pexeso_index->settings == NULL)
        return FAILED;

    pexeso_index->first_level = 0;
    pexeso_index->total_records = 0;

    pexeso_index->settings->root_directory = root_directory;
    pexeso_index->settings->vector_size = vector_size;
    pexeso_index->settings->max_dim_value = max_dim_value;
    pexeso_index->settings->min_dim_value = min_dim_value;
    pexeso_index->settings->leaf_cell_area = leaf_cell_area;

    return OK;
}

void exit_with_error(char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}