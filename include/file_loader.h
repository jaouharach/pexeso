#include <stdio.h>
#include <stdbool.h>

bool is_binaryfile(const char *filename);
response index_binary_files(pexeso_index * index, const char *bin_files_directory, unsigned int num_files, unsigned int base);
vector * load_binary_files(unsigned long * curr_total_vectors, const char *bin_files_directory, unsigned int num_files, unsigned int base, unsigned int num_dim);