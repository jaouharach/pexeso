#include<stdio.h>

struct vector{
	int table_id;
	int set_id;
	float * values;
};

struct cell {
    struct cell * parent;
    file_buffer * file_buffer;
    char * filename; 
    unsigned int is_leaf;
    vector * center;
    float edge_length;
    unsigned long num_vectors;
};


//In progress
void insert_vector(char * index_directory, vector *);
void append_vector_to_cell(cell * , vector *);

// done.
response init_leaf_cell(cell *, float length);
response init_leaf_cells(level * leaf_level, index_settings * settings);
/* initialize center vectors: ndc = number of distinct coordinates, k = vector length,  dim = vector length */
void init_center_vectors(float distinct_coordinates[], int ndc, int k, int dim, vector * center_vectors, vector temp, int append_at);



/* 
File name format : 
l[level]_center([center vector])_len[cell length]_size[num vectors].bin 
*/