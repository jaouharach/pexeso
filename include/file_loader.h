#include <stdio.h>
#include <stdbool.h>

bool is_binaryfile(const char *filename);
/* index raw binary vectors (in metric space) */
enum response index_binary_files(struct grid *grid, struct inv_index * index, const char *bin_files_directory, unsigned int num_files, unsigned int base);
/* index raw binary vectors (in metric space), only index a specific number af sets */
struct sid * index_query_binary_files(struct grid *grid, struct grid * Dgrid, struct inv_index * index, const char *bin_files_directory, unsigned int num_files, unsigned int base, int min_query_set_size, int max_query_set_size);
/* read all dataset files and create one big list of all the vectors in the dataset */
vector * load_binary_files(const char *bin_files_directory, unsigned long num_files, unsigned long long total_vectors, unsigned int base, unsigned int mtr_vector_length);

/* check datalake directory and count number of vectors (number of files is given as input */
unsigned long long get_datalake_info(const char *bin_files_directory, unsigned long num_files, unsigned long long *num_vectors, unsigned int *vector_length);
/* check datalake directory and count total files and number of vectors */
unsigned long long get_full_datalake_info(const char *bin_files_directory, unsigned long *num_files, unsigned long long *num_vectors, unsigned int *vector_length);
/* check query directory and count total files and number of vectors */
void get_query_data_info(const char *bin_files_directory, int num_query_sets, int min_query_set_size, int max_query_set_size, unsigned long *num_files, unsigned long long *num_vectors, unsigned int *vector_length);

/* save query results to csv file */
enum response  save_to_query_result_file(char * csv_file, struct sid * query_set, struct match_map * map);
/* make result file name and path */
char * make_file_path(char * work_dir, struct sid * query_set, unsigned int l, unsigned int dlsize, unsigned int vector_length, float runtime, unsigned int num_dist_calc, unsigned int total_checked_vec);
/* save query results to disk */
enum response save_results_to_disk(struct grid * Dgrid, struct grid * Qgrid, struct match_map * map);
/* create directory to store query results */
char * make_result_directory(char * work_dir, char* algorithm, unsigned int l, unsigned int num_query_sets, int min_query_set_size, int max_query_set_size);