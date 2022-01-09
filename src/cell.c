#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "../include/hgrid.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_buffer_manager.h"
#include "../include/file_buffer.h"
#include "../include/inv_index.h"

/* append vector to cell */
response append_vector_to_cell(struct grid *grid, struct inv_index * index, struct cell *cell,struct vector *vector, struct vector * v_mapping)
{
    // static int append_id = 0;
    // printf("append %d.\n\n\n", append_id+1);
    // append_id++;

    if (!get_file_buffer(grid, cell))
        exit_with_failure("Error in cell.c:  Could not get the \
                     file buffer for this cell.");

    // printf("Buffered list size = %u\n", cell->file_buffer->buffered_list_size);

    if (cell->file_buffer == NULL)
        exit_with_failure("Error in cell.c:  Couldn't append vector to cell NULL file buffer for \
                     this cell after creating it.");
    
    int mtr_vector_length = grid->settings->mtr_vector_length;
    int ps_vector_length = grid->settings->num_pivots;
    int max_leaf_size = grid->settings->max_leaf_size;
    int idx = cell->file_buffer->buffered_list_size;

    if (idx == 0) // if buffered list is empty
    {
        cell->file_buffer->mtr_buffered_list = NULL;
        cell->file_buffer->ps_buffered_list = NULL;
        
        cell->file_buffer->mtr_buffered_list = malloc(sizeof(struct v_type *) * max_leaf_size);
        cell->file_buffer->ps_buffered_list = malloc(sizeof(struct v_type *) * max_leaf_size);

        if (cell->file_buffer->mtr_buffered_list == NULL || cell->file_buffer->ps_buffered_list == NULL)
            exit_with_failure("Error in cell.c: Could not allocate memory for the buffered list.");
    }

    // get first address of memory block previously allocated for this cell
    cell->file_buffer->mtr_buffered_list[idx] = (v_type *) grid->buffer_manager->current_mtr_record;
    if (cell->file_buffer->mtr_buffered_list[idx] == NULL)
        exit_with_failure("Error in cell.c:  Could not \
                         allocate memory for the vector in the buffer.");

    // get first address of memory block previously allocated for this cell
    cell->file_buffer->ps_buffered_list[idx] = (v_type *) grid->buffer_manager->current_ps_record;
    if (cell->file_buffer->ps_buffered_list[idx] == NULL)
        exit_with_failure("Error in cell.c:  Could not \
                         allocate memory for the vector in the buffer.");


    // copy vector values to buffered list
    for (int i = 0; i < mtr_vector_length; ++i)
    {
        cell->file_buffer->mtr_buffered_list[idx][i] = vector->values[i];
    }
    for (int i = 0; i < ps_vector_length; ++i)
    {
        cell->file_buffer->ps_buffered_list[idx][i] = v_mapping->values[i];
    }

    grid->buffer_manager->current_mtr_record += sizeof(v_type) * mtr_vector_length;
    grid->buffer_manager->current_mtr_record_index++;

    grid->buffer_manager->current_mtr_record += sizeof(v_type) * ps_vector_length;
    grid->buffer_manager->current_ps_record_index++;

    cell->cell_size++;
    cell->file_buffer->buffered_list_size++;
    grid->total_records++;


    // track vector
    if (grid->settings->track_vector)
    {
        cell->vid[cell->cell_size - 1].table_id = vector->table_id;
        cell->vid[cell->cell_size - 1].set_pos = vector->set_id;
    }
    // add entry: set_id -> {cell} to inverted index    
    inv_index_append_entry(index, vector->table_id, vector->set_id, cell);

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
/* initialize cell */
response init_cell(cell *cell, float length, unsigned int num_child_cells)
{
    cell->edge_length = length;
    cell->num_child_cells = num_child_cells;
    cell->is_leaf = false;
    cell->cell_size = 0;
    
    // null pointers
    cell->parent = NULL;
    cell->children = NULL;
    cell->center = NULL;
    cell->file_buffer = NULL;
    cell->filename = NULL;
    cell->vid = NULL;

    return OK;
}

cell *get_child_cells(cell *parent_cell, unsigned int num_child_cells, level * children_level, struct grid_settings *settings)
{
    parent_cell->children = malloc(sizeof(struct cell) * num_child_cells);
    if (parent_cell->children == NULL)
        exit_with_failure("Error in cell.c: Couldn't allocate memory for child cells.");
    cell *child_cells = parent_cell->children;
    for (int c = 0; c < num_child_cells; c++)
    {
        // init child cell
        child_cells[c].parent = parent_cell;
        child_cells[c].num_child_cells = parent_cell->num_child_cells;
        child_cells[c].edge_length = parent_cell->edge_length / 2;
        child_cells[c].level = children_level;
        child_cells[c].center = NULL;
        child_cells[c].cell_size = 0;
        child_cells[c].file_buffer = NULL;
        child_cells[c].filename = NULL;
        child_cells[c].children = NULL;

        // printf("are leaf children = %s.\n", children_level->is_leaf ? "true" : "false");
        if(children_level->is_leaf)
        {
            child_cells[c].is_leaf = true;
            child_cells[c].vid = calloc(settings->max_leaf_size, sizeof(struct vid));
        }
        else
        {
            child_cells[c].is_leaf = false;
            child_cells[c].vid = NULL;
        }
        

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

    // free memory
    free(temp.values);
    int num_sign_vectors = (int) pow(2, settings->num_pivots);
    for(int s = num_sign_vectors - 1; s >= 0; s--)
        free(sign_arr[s].values);
    free(sign_arr);
    free(distinct_coor);
    
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
    dest->level = src->level;
    // printf("cell parent address = %p\n", dest->parent);
    dest->edge_length = src->edge_length;
    for (int i = 0; i < num_dim; i++)
        dest->center->values[i] = src->center->values[i];

    dest->is_leaf = src->is_leaf;
    dest->children = src->children;
    dest->num_child_cells = src->num_child_cells;
    dest->filename = src->filename;
    dest->cell_size = src->cell_size;
    dest->file_buffer = src->file_buffer;
    dest->vid = src->vid;

}
/* create cell filename */

enum response create_cell_filename(struct grid_settings *settings, struct cell * cell)
{
    /* 
        Each leaf cell has a file called: 
        level_edge_center(mag,mean)_numVectors
        02_0.334_(0.33,12.282)
        level: leaf level id (m)
        edge length: edge length of the leaf cell
        mag: magnitude of the center vector of the leaf cell
        num vectors: number of vectors stored in leaf cell
        number of punctuation marks (underscores:3, parentheses:2, commas:2): total 7
    */
    // if(cell->filename != NULL)
    // {
    //     printf("\n%s\n", cell->filename);
    //     exit_with_failure("Error in cell.c: Couldn't create cell filename, cell already has a filename!");
    // }

    int l = 0;
    cell->filename = malloc(sizeof(char) * (settings->max_filename_size));
    if(cell->filename == NULL)
        exit_with_failure("Error in cell.c: Couldn't allocate memory for cell filename.");

    l += sprintf(cell->filename+l ,"%02d", cell->level->id);
    // l += sprintf(cell->filename+l ,"%s", "_");
    l += sprintf(cell->filename+l ,"_%g", cell->edge_length);
    l += sprintf(cell->filename+l ,"_(%g,%g)", get_vector_magnitude(cell->center, settings->num_pivots), get_vector_mean(cell->center, settings->num_pivots));

    // printf("Cell filename = %s\n\n\n", cell->filename);
    return OK;
}

/* get list of vector in cell */
vector * get_vectors_mtr(struct cell * cell, unsigned int vector_length_mtr)
{
    if(!cell->is_leaf)
        exit_with_failure("Error in cell.c: Cannot get vectors of a non leaf cell!");

    if(cell->cell_size == 0)
        exit_with_failure("Error in cell.c: Cannot get vectors of an empty leaf cell!");

    // if file buffer not in disk return vectors in file buffer
    struct vector * cell_vectors = malloc(sizeof(struct vector) * cell->cell_size);
    for(int i = 0; i < cell->cell_size; i++)
        cell_vectors[i].values = malloc(sizeof(v_type) * vector_length_mtr);

    if(!cell->file_buffer->in_disk)
    {
        if(cell->cell_size != cell->file_buffer->buffered_list_size)
            exit_with_failure("Error in cell.c: Cell buffer not in disk yet cell_size != buffered list size!");
        
        // copy cell vectors into list of vectors
        for(int i = 0; i < cell->cell_size; i++)
        {
            for(int j = 0; j < vector_length_mtr; j++)
                cell_vectors[i].values[j] = cell->file_buffer->mtr_buffered_list[i][j];
                cell_vectors[i].set_id = cell->vid[i].set_pos;
                cell_vectors[i].table_id = cell->vid[i].table_id;
        }
    }
    // if file buffer is in disk load vectors from disk
    else
        exit_with_failure("Error in cell.c: can't get vectors, cell's file buffer is disk!");

    return cell_vectors;
}

/* get list of vector in cell (in pivot space) */
vector * get_vectors_ps(struct cell * cell, unsigned int vector_length_mtr)
{
    if(!cell->is_leaf)
        exit_with_failure("Error in cell.c: Cannot get vectors of a non leaf cell!");

    if(cell->cell_size == 0)
        exit_with_failure("Error in cell.c: Cannot get vectors of an empty leaf cell!");

    // if file buffer not in disk return vectors in file buffer
    struct vector * cell_vectors = malloc(sizeof(struct vector) * cell->cell_size);
    for(int i = 0; i < cell->cell_size; i++)
        cell_vectors[i].values = malloc(sizeof(v_type) * vector_length_mtr);

    if(!cell->file_buffer->in_disk)
    {
        if(cell->cell_size != cell->file_buffer->buffered_list_size)
            exit_with_failure("Error in cell.c: Cell buffer not in disk yet cell_size != buffered list size!");
        
        // copy cell vectors into list of vectors
        for(int i = 0; i < cell->cell_size; i++)
        {
            for(int j = 0; j < vector_length_mtr; j++)
                cell_vectors[i].values[j] = cell->file_buffer->ps_buffered_list[i][j];
                cell_vectors[i].set_id = cell->vid[i].set_pos;
                cell_vectors[i].table_id = cell->vid[i].table_id;
        }
    }
    // if file buffer is in disk load vectors from disk
    else
        exit_with_failure("Error in cell.c: can't get vectors, cell's file buffer is disk!");

    return cell_vectors;
}

/* get vector tuples (vectors in metric and ps space) */
// (v1, v1'), (v2, v2'), (v3, v3')
struct vector_tuple * get_vector_tuples(struct cell * cell, struct grid_settings * settings)
{
    if(!cell->is_leaf)
        exit_with_failure("Error in cell.c: Cannot get vectors of a non leaf cell!");

    if(cell->cell_size == 0)
        exit_with_failure("Error in cell.c: Cannot get vectors of an empty leaf cell!");

    // if file buffer not in disk return vectors in file buffer
    struct vector_tuple * cell_vectors = malloc(sizeof(struct vector_tuple) * cell->cell_size);
    for(int i = 0; i < cell->cell_size; i++)
    {
        cell_vectors[i].mtr_vector = malloc(sizeof(struct vector));
        cell_vectors[i].mtr_vector->values = malloc(sizeof(v_type) * settings->mtr_vector_length);
        cell_vectors[i].ps_vector->values = malloc(sizeof(v_type) * settings->num_pivots);

    }

    if(!cell->file_buffer->in_disk)
    {
        if(cell->cell_size != cell->file_buffer->buffered_list_size)
            exit_with_failure("Error in cell.c: Cell buffer not in disk yet cell_size != buffered list size!");
        
        // copy cell vectors into list of vector tuples
        for(int i = 0; i < cell->cell_size; i++)
        {
            for(int j = 0; j < settings->mtr_vector_length; j++)
                cell_vectors[i].mtr_vector->values[j] = cell->file_buffer->mtr_buffered_list[i][j];
            
            for(int j = 0; j < settings->num_pivots; j++)
                cell_vectors[i].mtr_vector->values[j] = cell->file_buffer->ps_buffered_list[i][j];
            
            cell_vectors[i].mtr_vector->set_id = cell->vid[i].set_pos;
            cell_vectors[i].mtr_vector->table_id = cell->vid[i].table_id;
            cell_vectors[i].ps_vector->set_id = cell->vid[i].set_pos;
            cell_vectors[i].ps_vector->table_id = cell->vid[i].table_id;
        }
    }
    // if file buffer is in disk load vectors from disk
    else
        exit_with_failure("Error in cell.c: can't get vectors, cell's file buffer is disk!");

    return cell_vectors;
}

/* get pointer to leaf cells of a given cell */
void get_leaf_cells(struct cell * cell, struct cell ** leaves, unsigned int * num_leaves)
{
    // if first call is with a leaf cell, leaf cell has no child cells
    if(cell->is_leaf && num_leaves == 0)
        exit_with_failure("Error in cell.c: Cannot find leaf cell for a leaf cell! a leaf cell has no child cells.");
    
    // a non leaf cell must have children
    if(cell->is_leaf == false && cell->num_child_cells == 0)
        exit_with_failure("Error in cell.c: Something went wrong, non leaf cell has no children.");

    // if a leaf cell is reached add it to list of leaves
    if(cell->is_leaf)
    {
        leaves = realloc(leaves, sizeof(struct cell *) * (*num_leaves + 1));
        leaves[*num_leaves] = cell;
        *num_leaves = *num_leaves + 1;

        return;
    }
    else
    {
        // go over all children of cell and collect lead cells 
        for(int i = 0; i < cell->num_child_cells; i++)
        {   
            get_leaf_cells(&cell->children[i], leaves, num_leaves);
        }
    }
}