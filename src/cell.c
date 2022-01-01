#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_buffer.h"

/* Find child cell with the closest center to vector */
cell *cell_route_to_closest_child(cell *parent_cell, vector *vector, unsigned int num_dim)
{
    float bsf = FLT_MAX;
    cell *closest_child_cell = NULL;
    for (int i = 0; i < parent_cell->num_child_cells; i++)
    {
        float d = euclidean_distance(vector, parent_cell->children[i].center, num_dim);
        if (d <= bsf)
        {
            bsf = d;
            closest_child_cell = &parent_cell->children[i];
        }
    }
    return closest_child_cell;
}

// /* initialize leaf cells in a level */
// response init_leaf_cells(level *leaf_level, index_settings *settings)
// {
//     /* allocate memory for cells */
//     leaf_level->cells = (cell *) malloc(sizeof(struct cell) * leaf_level->num_cells);
//     int i, j, ndc = log(leaf_level->num_cells) / log(settings->num_dim); //ndc = number of distinct coordinate values of cell center vectors

//     /* copmute distinct coordinate values (v1, v2, v3, ...) */
//     v_type * distinct_coordinates = (v_type *) malloc(sizeof(v_type) * ndc);
//     for (i = 0, j = 1; i < ndc; i++, j += 2)
//     {
//         distinct_coordinates[i] = settings->max_coordinate - ((j * settings->leaf_cell_edge_length) / 2);
//     }

//     /* init cells with their center vectors */
//     vector temp;
//     temp.values = malloc(sizeof(v_type) * settings->num_dim);

//     if (temp.values == NULL)
//         exit_with_failure("Error in cell.c: Could not allocate memory for temp vector.\n");

//     vector *center_vectors = malloc(sizeof(struct vector) * leaf_level->num_cells);

//     if (center_vectors == NULL)
//         exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");

//     for (i = 0; i < leaf_level->num_cells; i++)
//     {
//         center_vectors[i].values = malloc(sizeof(v_type) * settings->num_dim);

//         if (center_vectors[i].values == NULL)
//             exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");

//     }

//     init_center_vectors(distinct_coordinates, ndc, settings->num_dim,
//                         settings->num_dim, center_vectors, temp, 0);

//     for (i = 0; i < leaf_level->num_cells; i++)
//     {
//         init_leaf_cell(&leaf_level->cells[i], leaf_level->cell_edge_length);
//         leaf_level->cells[i].center = &center_vectors[i];

//         /* printf("(");
//         for(j = 0; j < settings->num_dim; j++)
//             printf("%.2f, ", leaf_level->cells[i].center->values[j]);
//         printf(")\n"); */
//     }

//     free(temp.values);
//     free(distinct_coordinates);
//     // don't free centre_vectors
//     return OK;
// }

// /* initialize cells in a level */
// response init_cells(level * level, index_settings *settings)
// {
//     /* allocate memory for cells */
//     level->cells = (cell *) malloc(sizeof(struct cell) * level->num_cells);
//     int i, j, ndc = log(level->num_cells) / log(settings->num_dim); //ndc = number of distinct coordinate values of cell center vectors

//     /* compute distinct coordinate values (v1, v2, v3, ...) */
//     v_type * distinct_coordinates = (v_type *) malloc(sizeof(v_type) * ndc);
//     for (i = 0, j = 1; i < ndc; i++, j += 2)
//     {
//         distinct_coordinates[i] = settings->max_coordinate - ((j * settings->leaf_cell_edge_length) / 2);
//     }

//     /* init cells with their center vectors */
//     vector temp;
//     temp.values = malloc(sizeof(v_type) * settings->num_dim);

//     if (temp.values == NULL)
//         exit_with_failure("Error in cell.c: Could not allocate memory for temp vector.\n");

//     vector *center_vectors = malloc(sizeof(struct vector) * level->num_cells);

//     if (center_vectors == NULL)
//         exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");

//     for (i = 0; i < level->num_cells; i++)
//     {
//         center_vectors[i].values = malloc(sizeof(v_type) * settings->num_dim);

//         if (center_vectors[i].values == NULL)
//             exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");

//     }

//     init_center_vectors(distinct_coordinates, ndc, settings->num_dim,
//                         settings->num_dim, center_vectors, temp, 0);

//     for (i = 0; i < level->num_cells; i++)
//     {
//         init_leaf_cell(&level->cells[i], level->cell_edge_length);
//         level->cells[i].center = &center_vectors[i];

//         /* printf("(");
//         for(j = 0; j < settings->num_dim; j++)
//             printf("%.2f, ", leaf_level->cells[i].center->values[j]);
//         printf(")\n"); */
//     }

//     free(temp.values);
//     free(distinct_coordinates);
//     // don't free centre_vectors
//     return OK;
// }

/* initialize leaf cell */
response init_leaf_cell(cell *cell, float length)
{
    cell->parent = NULL;
    cell->children = NULL;
    cell->num_child_cells = 0;

    cell->filename = malloc(sizeof(char));
    cell->file_buffer = NULL;

    cell->is_leaf = true;
    cell->edge_length = length;
    cell->center = NULL;

    file_buffer_init(cell);

    return OK;
}

/* initialize non leaf cell */
response init_cell(cell *cell, float length, unsigned int num_child_cells)
{
    cell->edge_length = length;
    cell->num_child_cells = num_child_cells;
    cell->is_leaf = false;

    // null pointers 
    cell->parent = NULL;
    cell->children = NULL;
    cell->center = NULL;
    cell->file_buffer = NULL;
    cell->filename = NULL;
    
    return OK;
}

cell * get_child_cells(cell * parent_cell, unsigned int num_child_cells, index_settings * settings)
{
    parent_cell->children = malloc(sizeof(struct cell) * num_child_cells);
    if(parent_cell->children == NULL)
        exit_with_failure("Error in cell.c: Couldn't allocate memory for child cells.");
    cell * child_cells = parent_cell->children;
    for(int c = 0; c < num_child_cells; c++)
    {
        child_cells[c].parent = parent_cell;
        child_cells[c].num_child_cells = parent_cell->num_child_cells;
        child_cells[c].edge_length = parent_cell->edge_length / 2;
        child_cells[c].is_leaf = false;
        child_cells[c].center = NULL;
    }
    // the values of center vectors for child cells will vary (from parent center vector) in earch direction d(., pi) by (+/-) parent_cell->edge_length /4
    unsigned int ndc = settings->num_dim * 2;
    v_type * distinct_coor = malloc(sizeof(v_type) * ndc); // 2 for (+) edge_length/4 and (-) edge_length/4
    if(distinct_coor == NULL)
        exit_with_failure("Error in cell.c: Couldn't allocate memory for distinct coordinate of child cells.");
    for(int i = 0, j = 0; i < settings->num_dim; i++, j += 2)
    {
        distinct_coor[j] =  parent_cell->center->values[i] + parent_cell->edge_length /4;
        distinct_coor[j+1] =  parent_cell->center->values[i] - parent_cell->edge_length /4;
    }

    // make center vectors for child cells.
    /* create center vectors */
    vector temp;
    temp.values = malloc(sizeof(v_type) * settings->num_dim);
    if (temp.values == NULL)
        exit_with_failure("Error in cell.c: Could not allocate memory for temp vector.\n");
    

    // allocate memory for center vectors
    vector *center_vectors = malloc(sizeof(struct vector) * parent_cell->num_child_cells);
    if (center_vectors == NULL)
        exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");

    for(int i = 0; i < parent_cell->num_child_cells; i++)
    {
        center_vectors[i].values = malloc(sizeof(v_type) * settings->num_dim);
        if (center_vectors[i].values == NULL)
            exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");
        
        parent_cell->children[i].center = &center_vectors[i];
    }


    // init_center_vectors(distinct_coor, ndc, settings->num_dim,
    //                     settings->num_dim, center_vectors, temp, 0);
    int s [] = {1, -1};
    vector * sign_arr = self_cartesian_product(s, settings->num_dim); // (1, 1, 1), (1, 1, -1), (1, -1, 1), ...
    for(int i = 0; i < parent_cell->num_child_cells; i++)
    {
        for(int j = 0; j < settings->num_dim; j++)
        {
            center_vectors[i].values[j] = parent_cell->center->values[j] + ((sign_arr[i].values[j] * parent_cell->edge_length) / 4);
        }
    }
    
    return parent_cell->children;
}
/* create center vectors for all cells of the data space */
void init_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector *center_vectors, vector temp, int append_at)
{
    static int curr_vector = 0;
    if (ndc == 0)
    {
        exit_with_failure("Error in cell.c: Couldn'l initialize center vector!\n");
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

/* initialize cells in a level */
// response make_level_cells(level *level, index_settings *settings)
// {
//     /* allocate memory for cells */
//     level->cells = (struct cell *)malloc(sizeof(struct cell) * level->num_cells);
//     if (level->cells == NULL)
//         exit_with_failure("Error in main.c: Could not allocate memory for next level cells.");
//     // initialize cells.
//     for (int c = 0; c < level->num_cells; c++)
//     {
//         init_cell(&level->cells[c], level->cell_edge_length);
//     }


//     /* find distinct center coordinates */
//     int i, j, ndc = settings->pivot_space_extrimity->values[0] / level->cell_edge_length; //ndc = number of distinct coordinate values of cell center vectors
//     printf("nuber of distinct coordinates = %d\n", ndc);
//     /* compute distinct coordinate values (v1, v2, v3, ...) */
//     v_type * distinct_coordinates = (v_type *) malloc(sizeof(v_type) * ndc);
//     for (i = 0, j = 1; i < ndc; i++, j += 2)
//     {
//         distinct_coordinates[i] = settings->pivot_space_extrimity->values[0] - ((j * level->cell_edge_length) / 2);
//     }

//     printf("Distinct coordinates:\n");
//     printf("(");
//     for(i = 0; i < ndc; i++)
//     {
//         printf("%f, ", distinct_coordinates[i]);
//     }
//     printf(")\n");


//     /* init cells with their center vectors */
//     vector temp;
//     temp.values = malloc(sizeof(v_type) * settings->num_dim);

//     if (temp.values == NULL)
//         exit_with_failure("Error in cell.c: Could not allocate memory for temp vector.\n");

//     vector *center_vectors = malloc(sizeof(struct vector) * level->num_cells);

//     if (center_vectors == NULL)
//         exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");

//     for (i = 0; i < level->num_cells; i++)
//     {
//         center_vectors[i].values = malloc(sizeof(v_type) * settings->num_dim);

//         if (center_vectors[i].values == NULL)
//             exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");
//     }

//     init_center_vectors(distinct_coordinates, ndc, settings->num_dim,
//                         settings->num_dim, center_vectors, temp, 0);

//     for (i = 0; i < level->num_cells; i++)
//     {
//         init_leaf_cell(&level->cells[i], level->cell_edge_length);
//         level->cells[i].center = &center_vectors[i];

//         /* printf("(");
//         for(j = 0; j < settings->num_dim; j++)
//             printf("%.2f, ", leaf_level->cells[i].center->values[j]);
//         printf(")\n"); */
//     }

//     free(temp.values);
//     free(distinct_coordinates);
//     // don't free centre_vectors
//     return OK;
// }
void cell_cpy(cell *dest, cell *src, unsigned int num_dim)
{
    dest->parent = src->parent;
    dest->edge_length = src->edge_length;
    for(int i = 0; i < num_dim; i++)
        dest->center->values[i] = src->center->values[i];

    dest->children = NULL;
    dest->num_child_cells = src->num_child_cells;
    dest->filename = src->filename;
    dest->file_buffer = src->file_buffer;

    dest->is_leaf = src->is_leaf;
}