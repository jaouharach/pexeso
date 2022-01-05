#include <stdio.h>
#include <float.h>
#include <dirent.h>
#include "globals.h"

struct index_settings {
    const char * root_directory;
    vector * pivots_mtr; // pivot vectors in metric space  
    vector * pivots_ps; // pivot vectors in pivot space  
    unsigned int num_pivots; // |P|, number of dimensions in pivot space
    vector * pivot_space_extremity;
    float pivot_space_volume;
    float leaf_cell_edge_length;
    unsigned int num_leaf_cells; // 2^(|P| * m)
    unsigned int num_levels; // m 
    unsigned int max_num_child_cells; // max number of child cells per parent cell
    unsigned int base; // number of bits per vector value in binary files (32, 64, ...)
    unsigned int mtr_vector_length; // vector length (number of dimensions) in metric space
    unsigned int vector_size; // vector size = (base/8) x vector length (in bytes)
    unsigned int max_filename_size; // each leaf cell is stored in one file 
    double buffered_memory_size;
    unsigned int max_leaf_size; // max number of vectors stored in one leaf cell
};

struct pexeso_index{
  unsigned long long total_records; // total vectors to be indexed
  struct level * first_level;
  struct index_settings * settings;
  struct file_buffer_manager * buffer_manager;
};

response init_index(const char * root_directory,
                unsigned int num_pivots,
                vector * extremity,
                unsigned int num_levels,
                unsigned long long num_vectors,
                unsigned int base,
                unsigned int mtr_vector_length,
                double buffered_memory_size,
                unsigned int max_leaf_size,
                pexeso_index * index);

vector * get_extremity(vector * pivot_vectors, unsigned int num_dim);

response index_insert(pexeso_index *, vector *);
void print_index(pexeso_index * index);