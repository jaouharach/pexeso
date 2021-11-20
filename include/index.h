#include<stdio.h>
#include <float.h>
#include"globals.h"

struct index_settings {
    const char * root_directory;  
    unsigned int num_dim;
    float max_coordinate;    
    float min_coordinate;
    float leaf_cell_edge_length;
    unsigned int num_levels;
};

struct index{
  unsigned long long total_records;
  struct level * first_level;
  struct index_settings * settings;
};

// done.
response init_index(const char * root_directory,
                unsigned int num_dim,
                float max_coordinate,    
                float min_coordinate,
                float leaf_cell_edge_length,
                index * index);

response append_vector(index * , vector *);
response index_binary_files(char * dataset_dir, unsigned int l, index *);