#include <stdio.h>

// run pexeso for one query column
enum response pexeso(struct sid * query_set, struct vector * query_vectors, unsigned int num_query_vectors,
                                struct grid * Dgrid, struct inv_index * inv_index);

/* quick browsing: evaluate leaf cells in Qgrid in Dgrid inverted index and get candidate pairs */
void quick_browse(struct grid * Dgrid, struct grid * Qgrid, struct pairs * pairs, struct grid_settings * settings);

/* create one query grid  */
enum response build_query_grid(struct grid * Qgrid, struct grid * Dgrid, inv_index * inv_index, struct vector * query_vectors, unsigned int num_query_vectors);