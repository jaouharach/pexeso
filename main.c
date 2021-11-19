#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "headers/main.h"

int main()
{
    const char *root_directory = "/home/jaouhara/Projects/pexeso/index";
    unsigned int num_dim = 3;
    float max_coordinate = 6; // we assume than all dims have the same min and max coordinate
    float min_coordinate = 0;
    float leaf_cell_edge_length = 2;
    int curr_state = OK;

    // initialize index
    index *pexeso_index = (index *)malloc(sizeof(struct index));

    if (pexeso_index == NULL)
        exit_with_error("Couldn't allocate memory for index!");

    if (!init_index(root_directory, num_dim, max_coordinate, min_coordinate, leaf_cell_edge_length, pexeso_index))
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

/* initialize leaf cells in a level */
int init_leaf_cells(level * leaf_level, index_settings * settings)
{
    /* allocate memory for cells */
    leaf_level->cells = (cell *)malloc(sizeof(struct cell) * leaf_level->num_cells);
    int i, j, ndc = log(leaf_level->num_cells)/log(settings->num_dim); //ndc = number of distinct coordinate values of cell center vectors
     
    /* copmute distinct coordinate values (v1, v2, v3, ...) */
    float *distinct_coordinates = (float *)malloc(sizeof(float) * ndc);
    for (i = 0, j = 1; i < ndc; i++, j += 2)
    {
        distinct_coordinates[i] = settings->max_coordinate - ((j * settings->leaf_cell_edge_length)/2);
    }

    /* init cells with their center vectors */
    vector temp;
    temp.values = malloc(sizeof(float) * settings->num_dim);

    vector * center_vectors  = malloc(sizeof(struct vector) * leaf_level->num_cells);
    for (i = 0; i < leaf_level->num_cells; i++)
    {
        center_vectors[i].values = malloc(sizeof(float) * settings->num_dim);
    }
    
    init_center_vectors(distinct_coordinates, ndc, settings->num_dim, 
                        settings->num_dim, center_vectors, temp, 0);

    for(i = 0; i < leaf_level->num_cells; i++)
    {
        init_leaf_cell(&leaf_level->cells[i], leaf_level->cell_length);
        leaf_level->cells[i].center = &center_vectors[i];

        printf("(");
        for(j = 0; j < settings->num_dim; j++)
            printf("%.2f, ", leaf_level->cells[i].center->values[j]);
        printf(")\n");
    }

    return OK;
}

/* initialize leaf cell */
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


/* print error and kill process */
void exit_with_error(char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

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


/* initialize leaf level */
int init_leaf_level(index_settings *settings, level *level)
{
    //num_cells = (max - min) ^ num_dim /cell_length ^ num_dim
    
    // check if num cells is integer
    if (fmod(
            pow(settings->max_coordinate - settings->min_coordinate, settings->num_dim),
            pow(settings->leaf_cell_edge_length, settings->num_dim) == 0))
    {
        level->id = 0;
        level->num_cells = (unsigned int)abs(
            pow(settings->max_coordinate - settings->min_coordinate, settings->num_dim) / pow(settings->leaf_cell_edge_length, settings->num_dim));
        
        printf("num cells: %d\n", level->num_cells);

        level->cell_length = settings->leaf_cell_edge_length;
        level->next_level = NULL;

        // initialize cells.
        init_leaf_cells(level, settings);

        return OK;
    }
    else
        exit_with_error("Warning: Number of cells can only be integer! please change settings.");
}

void init_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector * center_vectors, vector temp, int append_at)
{
    static int curr_vector = 0;
    if (ndc == 0 || k > ndc)
    {
        exit_with_error("Couldn'l initialize center vector!\n");
    }

    if (k == 0)
    {
        // add temp to sub_array
        for(int x = 0; x < dim; x++)
            center_vectors[curr_vector].values[x] = temp.values[x];
        curr_vector++;
        return;
    }

    for (int j = 0; j < ndc; j++)
    {
        // add distinct_coordinates[j] to temp values
        temp.values[append_at] = distinct_coordinates[j];
        init_center_vectors(distinct_coordinates, ndc, k - 1, dim, center_vectors, temp, append_at + 1);;
    }
}