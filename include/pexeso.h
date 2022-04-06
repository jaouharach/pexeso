#include <stdio.h>


/* pexeso set similarity search algorithm */
void pexeso(const char * query_file_dir, struct grid * Dgrid, struct inv_index * inv_index);

/* create query grid */
struct sid * build_query_grid(struct grid * Qgrid, struct grid * Dgrid, inv_index * inv_index, const char * query_file_dir, long long unsigned int * num_query_vectors);