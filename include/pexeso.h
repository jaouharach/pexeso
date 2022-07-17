#include <stdio.h>


/* pexeso set similarity search algorithm */
void pexeso(const char * query_file_dir, struct grid * Dgrid, struct inv_index * inv_index);

/* quick browsing: evaluate leaf cells in Qgrid in Dgrid inverted index and get candidate pairs */
void quick_browse(struct grid * Dgrid, struct grid * Qgrid, struct pairs * pairs, struct grid_settings * settings);

/* create query grid */
struct sid * build_query_grid(struct grid * Qgrid, struct grid * Dgrid, inv_index * inv_index, const char * query_file_dir, long long unsigned int * num_query_vectors);