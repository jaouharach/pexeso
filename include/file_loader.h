#include <stdio.h>
#include <stdbool.h>

bool is_binaryfile(const char *filename);
/* index raw binary vectors (in metric space) */
response index_binary_files(pexeso_index *index, const char *bin_files_directory, unsigned int num_files, unsigned int base);
/* read all dataset files and create one big list of all the vectors in the dataset */
vector * load_binary_files(const char *bin_files_directory, unsigned long num_files, unsigned long long total_vectors, unsigned int base, unsigned int mtr_vector_length);
/* check dataset directory and count total files and number of vectors */
unsigned long long get_dataset_info(const char *bin_files_directory, unsigned long *num_files, unsigned long long *num_vectors, unsigned int *vector_length);
