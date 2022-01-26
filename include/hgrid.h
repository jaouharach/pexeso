#include <stdio.h>
#include <float.h>
#include <dirent.h>
#include "globals.h"

struct grid_settings {
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
    double mtr_buffered_memory_size; // memory array size for metric space memory array
    double ps_buffered_memory_size; // memory array size for pivot space memory array
    unsigned int max_leaf_size; // max number of vectors stored in one leaf cell
    unsigned int track_vector; // (value = 0 / 1) track data vector id or not, vector id = (table_id, column_id) 
    struct query_settings * query_settings; // distance threshold tau, joinability threshold T ...
};

/* hierarchical grid structure */
struct grid{
  unsigned long long total_records; // total vectors to be indexed
  struct level * root; // root level, num cells  = 1.
  struct level * first_level;
  struct grid_settings * settings;
  struct file_buffer_manager * buffer_manager;
  bool is_query_grid;// if the grid is a Qgrid
  //file to track vector ids 
  const char * vid_filename; 
  FILE * vid_file;
  unsigned int vid_pos_ctr;
  struct vid * vid_cache;
};

enum response init_grid(const char *root_directory,
                    unsigned int num_pivots,
                    vector *pivots_mtr,
                    vector *pivots_ps,
                    vector *extremity,
                    unsigned int num_levels,
                    unsigned long long num_vectors,
                    unsigned int base,
                    unsigned int mtr_vector_length,
                    double mtr_buffered_memory_size,
                    double ps_buffered_memory_size,
                    unsigned int max_leaf_size,
                    unsigned int track_vector,
                    bool is_query_grid,
                    struct query_settings * query_settings,
                    struct grid *grid);

/* get the farthest vector in the pivot space */
vector * get_extremity(vector * pivot_vectors, unsigned int num_dim);

/* insert a vector in the grid */
enum response grid_insert(struct grid *, struct inv_index *, vector *);

/* print grid in console */
void dump_grid_to_console(struct grid * grid);

/* write grid to disk */
enum response grid_write(struct grid *grid);

/* destroy grid */
enum response grid_destroy(struct grid *grid);

/* destroy grid levels*/
enum response grid_destroy_level(struct grid *grid, struct level *level);

/* destroy buffer manager */
enum response destroy_buffer_manager(struct grid *grid);

/* write level to disk */
enum response level_write(struct grid *grid, struct level *level, FILE *file);

/* destroy query grid */
enum response query_grid_destroy(struct grid *grid);