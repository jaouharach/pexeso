#include <stdio.h>
#include <float.h>
#include <dirent.h>
#include "globals.h"

struct index_settings {
    const char * root_directory;  
    unsigned int num_dim; // P
    float max_coordinate;    
    float min_coordinate;
    float leaf_cell_edge_length;
    unsigned int num_leaf_cells; // 2^(P * m)
    unsigned int num_levels; // m 
};

struct pexeso_index{
  unsigned long long total_records;
  struct level * first_level;
  struct index_settings * settings;
};

// done.
response init_index(const char * root_directory,
                unsigned int num_dim,
                float max_coordinate,    
                float min_coordinate,
                unsigned int num_levels,
                pexeso_index * index);

response insert_vector(pexeso_index *, vector *);
void display_indexed_vectors(pexeso_index *);