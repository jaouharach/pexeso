#include<stdio.h>
#include"globals.h"

struct index_settings {
    const char * root_directory;  
    unsigned int num_dim;
    float max_coordinate;    
    float min_coordinate;
    float leaf_cell_edge_length;
};

struct index{
  unsigned long long total_records;
  level * first_level;
  struct index_settings * settings;
};

//In progress
void create_index(char * index_directory, char * dataset_directory);

// done.
int init_index(const char * root_directory,
                unsigned int num_dim,
                float max_coordinate,    
                float min_coordinate,
                float leaf_cell_edge_length,
                index * pexeso_index);

// file for each cell