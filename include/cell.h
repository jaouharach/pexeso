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
void create_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector * center_vectors, vector temp, int append_at);
cell * cell_route_to_closest_child (cell * parent_cell, vector * vector, unsigned int num_dim);
response init_cell(cell *cell, float length, unsigned int num_child_cells); // initialize non leaf cell
/* create child cells for a parent cell */
cell * get_child_cells(cell *parent_cell, unsigned int num_child_cells, level * children_level, struct grid_settings *settings);
void cell_cpy(cell *dest, cell *src, unsigned int num_dim);
/* append vector to cell */
response append_vector_to_cell(struct grid *grid, struct inv_index * index, struct cell *cell,struct vector *vector);
/* create cell filename */
enum response create_cell_filename(struct grid_settings *settings, struct cell * cell);
/* get list of vectors stored in cell */
vector * get_vectors(struct cell * cell, unsigned int num_pivots);
/* get pointer to leaf cells of a given cell */
void get_leaf_cells(struct cell * cell, struct cell ** leaves, unsigned int * num_leaves);