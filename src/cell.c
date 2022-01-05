#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_buffer_manager.h"
#include "../include/file_buffer.h"

/* append vector to cell */
response append_vector_to_cell(struct pexeso_index *index, struct cell *cell,struct vector *vector)
{
    if (!get_file_buffer(index, cell))
        exit_with_failure("Error in cell.c:  Could not get the \
                     file buffer for this cell.");

    
    if (cell->file_buffer == NULL)
        exit_with_failure("Error in cell.c:  Couldn't append vector to cell NULL file buffer for \
                     this cell after creating it.");
    
    int idx = cell->file_buffer->buffered_list_size;
    int vector_length = index->settings->mtr_vector_length;
    int max_leaf_size = index->settings->max_leaf_size;

    if (idx == 0) // if buffered list is empty
    {
        cell->file_buffer->buffered_list = NULL;
        cell->file_buffer->buffered_list = malloc(sizeof(struct ts_type *) * max_leaf_size);

        if (cell->file_buffer->buffered_list == NULL)
            exit_with_failure("Error in cell.c:  Could not \
                                allocate memory for the buffered list.");
    }

    // get first address of memory block previously allocated for this cell
    cell->file_buffer->buffered_list[idx] = (v_type *) index->buffer_manager->current_record;
    index->buffer_manager->current_record += sizeof(v_type) * vector_length;
    index->buffer_manager->current_record_index++;

    if (cell->file_buffer->buffered_list[idx] == NULL)
        exit_with_failure("Error in cell.c:  Could not \
                         allocate memory for the vector in the buffer.");

    // copy vector values to buffered list
    for (int i = 0; i < vector_length; ++i)
    {
        cell->file_buffer->buffered_list[idx][i] = vector->values[i];
    }

    // (todo) track vector
    // cell->vid[cell->cell_size - 1].table_id = vector->table_id;
    // cell->vid[cell->cell_size - 1].set_id = vector->set_id;


    ++cell->file_buffer->buffered_list_size;
    index->buffer_manager->current_values_count += vector_length;

    return OK;
}
/* Find child cell with the closest center to vector */
cell *cell_route_to_closest_child(cell *parent_cell, vector *vector, unsigned int num_dim)
{
    float bsf = FLT_MAX;
    cell *closest_child_cell = NULL;
    for (int i = 0; i < parent_cell->num_child_cells; i++)
    {
        float d = euclidean_distance(vector, parent_cell->children[i].center, num_dim);
        if (d < bsf)
        {
            bsf = d;
            closest_child_cell = &parent_cell->children[i];
        }
    }
    return closest_child_cell;
}

/* initialize leaf cell */
response to_leaf_cell(struct cell *cell, struct cell *parent, char *filename)
{
    cell->parent = parent;
    cell->children = NULL;
    cell->is_leaf = true;

    cell->filename = filename;

    cell->file_buffer = NULL;
    // file_buffer_init(cell);

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

cell *get_child_cells(cell *parent_cell, unsigned int num_child_cells, index_settings *settings)
{
    parent_cell->children = malloc(sizeof(struct cell) * num_child_cells);
    if (parent_cell->children == NULL)
        exit_with_failure("Error in cell.c: Couldn't allocate memory for child cells.");
    cell *child_cells = parent_cell->children;
    for (int c = 0; c < num_child_cells; c++)
    {
        child_cells[c].parent = parent_cell;
        child_cells[c].num_child_cells = parent_cell->num_child_cells;
        child_cells[c].edge_length = parent_cell->edge_length / 2;
        child_cells[c].is_leaf = false;
        child_cells[c].center = NULL;
    }
    // the values of center vectors for child cells will vary (from parent center vector) in earch direction d(., pi) by (+/-) parent_cell->edge_length /4
    unsigned int ndc = settings->num_pivots * 2;
    v_type *distinct_coor = malloc(sizeof(v_type) * ndc); // 2 for (+) edge_length/4 and (-) edge_length/4
    if (distinct_coor == NULL)
        exit_with_failure("Error in cell.c: Couldn't allocate memory for distinct coordinate of child cells.");
    for (int i = 0, j = 0; i < settings->num_pivots; i++, j += 2)
    {
        distinct_coor[j] = parent_cell->center->values[i] + parent_cell->edge_length / 4;
        distinct_coor[j + 1] = parent_cell->center->values[i] - parent_cell->edge_length / 4;
    }

    /* create center vectors for child cells */
    vector temp;
    temp.values = malloc(sizeof(v_type) * settings->num_pivots);
    if (temp.values == NULL)
        exit_with_failure("Error in cell.c: Could not allocate memory for temp vector.\n");

    // allocate memory for center vectors
    vector *center_vectors = malloc(sizeof(struct vector) * parent_cell->num_child_cells);
    if (center_vectors == NULL)
        exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");

    for (int i = 0; i < parent_cell->num_child_cells; i++)
    {
        center_vectors[i].values = malloc(sizeof(v_type) * settings->num_pivots);
        if (center_vectors[i].values == NULL)
            exit_with_failure("Error in cell.c: Could not allocate memory for list of center vectors.\n");

        parent_cell->children[i].center = &center_vectors[i];
    }

    /* 
        given parent_cell with center vector c = [ p1, p2, p3, ...] 
        create child center vectors 
        ch [p1 (+/-) edge_len/4,  p2 (+/-) edge_len/4, p3 (+/-) edge_len/4, ...]

    */

    int s[] = {1, -1};                                                  // sign (+/-)
    vector *sign_arr = self_cartesian_product(s, settings->num_pivots); // (1, 1, 1), (1, 1, -1), (1, -1, 1), ...
    for (int i = 0; i < parent_cell->num_child_cells; i++)
    {
        for (int j = 0; j < settings->num_pivots; j++)
        {
            center_vectors[i].values[j] = parent_cell->center->values[j] + ((sign_arr[i].values[j] * parent_cell->edge_length) / 4);
        }
    }

    return parent_cell->children;
}
/* create center vectors for all cells of the data space */
void create_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector *center_vectors, vector temp, int append_at)
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
        create_center_vectors(distinct_coordinates, ndc, k - 1, dim, center_vectors, temp, append_at + 1);
    }
}

void cell_cpy(cell *dest, cell *src, unsigned int num_dim)
{
    dest->parent = src->parent;
    dest->edge_length = src->edge_length;
    for (int i = 0; i < num_dim; i++)
        dest->center->values[i] = src->center->values[i];

    dest->children = NULL;
    dest->num_child_cells = src->num_child_cells;
    dest->filename = src->filename;
    dest->file_buffer = src->file_buffer;
}