#include <stdio.h>
#include <stdlib.h>

/* pivot selection algorithm */
vector * select_pivots(vector * dataset, int * dataset_dim, unsigned int num_pivots, unsigned int fft_scale);

/* FFT: Farthest First Traversal */
vector *fft(vector *data_set, int * dataset_dim, unsigned int k);
// min distance to a set of outliers (for fft)
float min_distance(vector *outliers, unsigned int num_outliers, vector *v, unsigned int v_len);
