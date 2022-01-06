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
    unsigned int track_vector; // (value = 0 / 1) track data vector id or not, vector id = (table_id, column_id) 
};

struct pexeso_index{
  unsigned long long total_records; // total vectors to be indexed
  struct level * first_level;
  struct index_settings * settings;
  struct file_buffer_manager * buffer_manager;

  //file to track vector ids 
  const char * vid_filename; 
  FILE * vid_file;
  unsigned int vid_pos_ctr;
  struct vid * vid_cache;
};

response init_index(const char * root_directory,
                unsigned int num_pivots,
                vector * pivots_mtr, 
                vector * pivots_ps,
                vector * extremity,
                unsigned int num_levels,
                unsigned long long num_vectors,
                unsigned int base,
                unsigned int mtr_vector_length,
                double buffered_memory_size,
                unsigned int max_leaf_size,
                unsigned int track_vector,
                pexeso_index * index);

/* get the farthest vector in the pivot space */
vector * get_extremity(vector * pivot_vectors, unsigned int num_dim);

/* insert a vector in the index */
response index_insert(pexeso_index *, vector *);

/* print index in console */
void print_index(pexeso_index * index);

/* write index to disk */
enum response index_write(pexeso_index *index);

/* destroy index */
void index_destroy(struct pexeso_index *index, struct level *level);

/* destroy buffer manager */
enum response destroy_buffer_manager(struct pexeso_index *index);