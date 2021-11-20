#include<stdio.h>

struct vector{
	int table_id;
	int set_id;
	float * values;
};

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


//In progress
void insert_vector(char * index_directory, vector *);

// done.
response init_leaf_cell(cell *, float length);
response init_leaf_cells(level * leaf_level, index_settings * settings);
/* initialize center vectors: ndc = number of distinct coordinates, k = vector length,  dim = vector length */
void init_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector * center_vectors, vector temp, int append_at);
float euclidean_distance(vector *, vector *, int k);
cell * cell_route_to_closest_child (cell * parent_cell, vector * vector, unsigned int num_dim);


/* 
File name format : 
l[level]_center([center vector])_len[cell length]_size[num vectors].bin 
*/