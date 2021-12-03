#include <stdio.h>
#include <float.h>
#include <dirent.h>
#include "globals.h"

struct index_settings {
    const char * root_directory;  
    unsigned int num_dim; // |P|
    vector * pivot_space_extrimity; // 
    float pivot_space_volume;
    float leaf_cell_edge_length;
    unsigned int num_leaf_cells; // 2^(|P| * m)
    unsigned int num_levels; // m 
};

struct pexeso_index{
  unsigned long long total_records;
  struct level * first_level;
  struct index_settings * settings;
};

response init_index(const char * root_directory,
                unsigned int num_pivots,
                vector * extrimity,
                unsigned int num_levels,
                pexeso_index * index);

vector * get_extrimity(vector * pivot_vectors, unsigned int num_dim);

response insert_vector(pexeso_index *, vector *);
void display_indexed_vectors(pexeso_index *);