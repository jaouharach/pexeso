#include <stdio.h>
#include <stdlib.h>

/* pivot selection algorithm */
vector * select_pivots(vector * dataset, int * dataset_dim, unsigned int num_pivots, unsigned int fft_scale);

/* FFT: Farthest First Traversal */
vector *fft(vector *data_set, int * dataset_dim, unsigned int k);
// min distance to a set of outliers (for fft)
float min_distance(vector *outliers, unsigned int num_outliers, vector *v, unsigned int v_len);

/* map vector to pivot space  v -> v_mapping*/
void map_vector(vector *v, unsigned int v_len, vector *v_mapping, vector *pivots, unsigned int num_pivots);

/*  get distance matrix of a dataset (mapping of vectors to the pivot space), 
    dataset_mtr (dataset in metric space) --> dataset_ps (dataset in pivot space) */
vector * map_to_pivot_space(vector * dataset_mtr, int * dataset_dim, vector * pivots, unsigned int num_pivots);

/* EMPCA: Expectation Management for Principal Component Analysis */
gsl_matrix * empca(vector *data_set, unsigned int num_vectors, unsigned int dim, unsigned int num_pc);
