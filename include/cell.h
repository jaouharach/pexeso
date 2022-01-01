#include <stdio.h>
#include <stdbool.h>


struct cell {
    struct cell * parent;
    struct cell * children;
    unsigned int num_child_cells;

    file_buffer * file_buffer;
    char * filename; 

    bool is_leaf;
    vector * center;
    float edge_length;
};


response to_leaf_cell(struct cell *cell, struct cell *parent, char * filename);
/* initialize center vectors: ndc = number of distinct coordinates, k = vector length,  dim = vector length */
void create_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector * center_vectors, vector temp, int append_at);
cell * cell_route_to_closest_child (cell * parent_cell, vector * vector, unsigned int num_dim);
response init_cell(cell *cell, float length, unsigned int num_child_cells); // initialize non leaf cell
/* create child cells for a parent cell */
cell * get_child_cells(cell * parent_cell, unsigned int num_child_cells, index_settings *settings);
void cell_cpy(cell *dest, cell *src, unsigned int num_dim);