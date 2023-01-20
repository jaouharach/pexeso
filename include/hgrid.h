#include <stdio.h>
#include <float.h>
#include <dirent.h>
#include "globals.h"
#include <string.h>
struct grid_settings {
    char * root_directory; // where cell files will be stored
    char * query_root_directory; // where query grid cell files will be stored
    char * work_directory; // where Hgrid and Qrid will be stored
    unsigned int fft_scale;
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

struct query_settings
{
    float join_threshold; // T 
    float dist_threshold; // tau
    int num_query_sets;
    int min_query_set_size;
    int max_query_set_size;
    unsigned int mtr_buffered_memory_size; // amount of memory (in MB) to store metric space vectors for query grid
    unsigned int ps_buffered_memory_size; // amount of memory (in MB) to store pivot space vectors for query grid
};

/* hierarchical grid structure */
struct grid{
  unsigned long long total_records; // total vectors to be indexed
  struct level * root; // root level, num cells  = 1.
  struct level * first_level;
  struct grid_settings * settings;
  struct stats_info * stats;
  struct file_buffer_manager * buffer_manager;
  bool is_query_grid;// if the grid is a Qgrid
  //file to track vector ids 
  const char * vid_filename; 
  FILE * vid_file;
  unsigned int vid_pos_ctr;
  struct vid * vid_cache;
};

struct stats_info {
  // counters
  unsigned long total_cells_count;
  unsigned long leaf_cells_count;
  unsigned int empty_leaf_cells_count;

  double loaded_files_size; // total size of loaded data tables (in Bytes)
  double loaded_query_files_size; // total size of loaded query tables (in Bytes)

  unsigned long loaded_files_count; // total loaded data tables 
  unsigned long loaded_query_files_count; // total loaded query tables 
  unsigned long loaded_sets_count; // total data sets (columns)
  unsigned long loaded_query_sets_count; // total query sets (columns)

  unsigned long loaded_vec_count; // total loaded data vectors 
  unsigned long loaded_qvec_count; // total loaded query vectors 


  unsigned long out_of_ps_space_vec_count; // data vectors that are out of the pivot space
  unsigned long out_of_ps_space_qvec_count; // query vectors that are out of the pivot space

  unsigned long total_queries_count;

  unsigned int used_lemmas_count[7]; // count how many times each lemma has been used

  unsigned int filtered_cells_count;
  unsigned int visited_cells_count;
  unsigned int visited_matching_cells_count;
  unsigned int visited_candidate_cells_count;

  unsigned long filtered_vectors_count;
  unsigned long matching_vectors_count;
  unsigned long checked_vectors_in_ps_count; // nb of vectors checked in the pivots space
  unsigned long checked_vectors_in_mtr_count; // nb of vectors checked in the metric space

  unsigned int count_add_mpair;
  unsigned int count_add_cpair;

  unsigned long count_dist_calc;
  
  // timers 
  double idx_append_vec_to_leaf_total_time;	
	double idx_append_vec_to_leaf_input_time;
	double idx_append_vec_to_leaf_output_time;
	double idx_append_vec_to_leaf_cpu_time;

  double total_input_time;
  double total_output_time;
  double total_parse_time;
  double total_query_time;
  double total_pivot_selection_time;
  double total_time;

  double grid_building_total_time; // time to build the grid from scratch
  double grid_building_input_time;
  double grid_building_output_time;
  double grid_writing_total_time; // time to write the grid to disk
};

enum response init_grid(const char *work_dir,
                    unsigned int fft_scale,
                    unsigned int num_pivots,
                    vector *pivots_mtr,
                    vector *pivots_ps,
                    vector *extremity,
                    unsigned int num_levels,
                    unsigned long long num_vectors,
                    unsigned int base,
                    unsigned int mtr_vector_length,
                    double mtr_buffered_memory_size,
                    unsigned int max_leaf_size,
                    unsigned int track_vector,
                    bool is_query_grid,
                    struct query_settings * query_settings,
                    struct grid *grid,
                    unsigned int read_from_disk);

/* init statistics */
enum response init_grid_stats(struct grid * grid);

// rest all query related stats
enum response reset_grid_stats(struct grid * grid);

/* get the farthest vector in the pivot space */
vector * get_extremity(vector * pivot_vectors, unsigned int num_dim);

/* insert a vector in the grid */
enum response grid_insert(struct grid *, struct inv_index *, vector *);

/* print grid in console */
void dump_grid_to_console(struct grid * grid);

/* write grid to disk */
enum response grid_write(struct grid *grid);

/* read grid from disk */
struct grid * grid_read(const char * grid_directory, struct query_settings * query_settings);

/* read level from disk */
enum response read_leaf_level(struct grid *grid, FILE *root_file, struct level * leaf_level);

/* destroy grid */
enum response grid_destroy(struct grid *grid);

/* destroy grid levels*/
enum response grid_destroy_level(struct grid *grid, struct level *level);

/* destroy buffer manager */
enum response destroy_buffer_manager(struct grid *grid);

/* write level to disk */
enum response write_leaf_level(struct grid *grid, struct level *level, FILE *file);

/* destroy query grid */
enum response query_grid_destroy(struct grid *grid);

/* print stats */
void print_grid_stats(struct grid * grid);

/* print settings */
void print_grid_settings(struct grid * grid);

/* count empty leaf cells */
enum response count_empty_leaf_cells(struct grid * grid);