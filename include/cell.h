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


response init_leaf_cell(cell *, float length);
response init_leaf_cells(level * leaf_level, index_settings * settings);
/* initialize center vectors: ndc = number of distinct coordinates, k = vector length,  dim = vector length */
void init_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector * center_vectors, vector temp, int append_at);
cell * cell_route_to_closest_child (cell * parent_cell, vector * vector, unsigned int num_dim);
response init_cell(cell *cell, float length); // initialize non leaf cell
/* create cells for a given level */
response make_level_cells(level * level, index_settings *settings);
