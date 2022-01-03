#include <stdio.h>
#include <float.h>
#include <dirent.h>
#include "globals.h"

struct index_settings {
    const char * root_directory;  
    unsigned int num_pivots; // |P|, number of dimensions in pivot space
    vector * pivot_space_extremity; // 
    float pivot_space_volume;
    float leaf_cell_edge_length;
    unsigned int num_leaf_cells; // 2^(|P| * m)
    unsigned int num_levels; // m 
    unsigned int max_num_child_cells; // max number of child cells per parent cell
    unsigned int base; // number of bits per vector value in binary files (32, 64, ...)
    unsigned int mtr_vector_length; // vector length (number of dimensions) in metric space
    unsigned int vector_size; // vector size = (base/8) x vector length (in bytes)
    unsigned int max_filename_size; // each leaf cell is stored in one file 
};

struct pexeso_index{
  unsigned long long total_records; // total vectors to be indexed
  struct level * first_level;
  struct index_settings * settings;
};

response init_index(const char * root_directory,
                unsigned int num_pivots,
                vector * extremity,
                unsigned int num_levels,
                unsigned long long num_vectors,
                unsigned int base,
                unsigned int mtr_vector_length,
                pexeso_index * index);

vector * get_extremity(vector * pivot_vectors, unsigned int num_dim);

response insert_vector(pexeso_index *, vector *);
void display_indexed_vectors(pexeso_index *);