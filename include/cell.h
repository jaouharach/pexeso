#include <stdio.h>
#include <stdbool.h>


struct cell {
    struct level * level;
    struct cell * parent;
    struct cell * children;
    unsigned int num_child_cells;

    struct file_buffer * file_buffer;
    char * filename; 

    bool is_leaf;
    unsigned int cell_size; // number of vectors currently stored in cell
    
    vector * center;
    float edge_length;

    unsigned int vid_pos; // first position in vid.idx file
    struct vid * vid;

    
};

/* initialize center vectors: ndc = number of distinct coordinates, k = vector length,  dim = vector length */
void create_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector * center_vectors, vector temp, int append_at, int *curr_vector);

cell * cell_route_to_closest_child (cell * parent_cell, vector * vector, unsigned int num_dim);

response init_cell(cell *cell, float length, unsigned int num_child_cells); // initialize non leaf cell

/* create child cells for a parent cell */
cell * get_child_cells(cell *parent_cell, unsigned int num_child_cells, level * children_level, struct grid_settings *settings);
void cell_cpy(cell *dest, cell *src, unsigned int num_dim);

/* append vector to cell */
response append_vector_to_cell(struct grid *grid, struct inv_index * index, struct cell *cell,struct vector *vector, struct vector * v_mapping);

/* create cell filename */
enum response create_cell_filename(struct grid_settings *settings, struct cell * cell);

/* get list of vector in cell */
vector * get_vectors_mtr(struct cell * cell, unsigned int vector_length_mtr);

/* get list of vector in cell (in pivot space) */
vector * get_vectors_ps(struct cell * cell, unsigned int vector_length_mtr);

/* get pointer to leaf cells of a given cell */
void get_leaf_cells(struct cell * cell, struct cell ** leaves, unsigned int * num_leaves);

/* get number of leaf cells of a given cell */
void get_num_leaf_cells(struct cell * cell, unsigned int * num_leaves);

/* get vector tuples (vectors in metric and ps space) */
// (v1, v1'), (v2, v2'), (v3, v3')
struct vector_tuple * get_vector_tuples(struct cell * cell, struct grid_settings * settings);

/* get list of vector in the sub leaf cells of a non leaf cell (in pivot space) */
vector * get_sub_cells_vectors_ps(struct cell * cell, unsigned int num_pivots, long unsigned int * num_vectors);