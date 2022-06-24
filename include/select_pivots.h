#include <stdio.h>
#include <stdlib.h>

struct best_fft
{
    int fft_scale;
    float exremity;
    vector * pivots_mtr;
};

/* pivot selection algorithm */
vector * select_pivots(vector * dataset, int * dataset_dim, unsigned int num_pivots, unsigned int fft_scale);

/* FFT: Farthest First Traversal */
vector *fft(vector *data_set, int * dataset_dim, unsigned int k);
// min distance to a set of outliers (for fft)
float min_distance(vector *outliers, unsigned int num_outliers, vector *v, unsigned int v_len);

/* map vector to pivot space  v -> v_mapping*/
enum response map_vector(struct vector *v, unsigned int v_len, struct vector *v_mapping, struct grid_settings * settings, bool in_query_grid);

/*  get distance matrix of a dataset (mapping of vectors to the pivot space), 
    dataset_mtr (dataset in metric space) --> dataset_ps (dataset in pivot space) */
vector * map_to_pivot_space(vector * dataset_mtr, int * dataset_dim, vector * pivots, unsigned int num_pivots);

/* EMPCA: Expectation Management for Principal Component Analysis */
gsl_matrix * empca(vector *data_set, unsigned int num_vectors, unsigned int dim, unsigned int num_pc);

/* select pivots from the empca result (pc with highest projection on dataset)
    returns indecies of the pivots with the highest projection */
int * select_pivots_by_pca_result_angle(gsl_matrix * pca_result, int *result_dim, int num_pivots);

/* check if a new outlier is already in the list of outliers */
unsigned int in_outliers(int * outliers_idx, int new_outlier, int num_outliers);

/* select best fft scale */
struct vector *select_pivots_with_best_fft_scale(struct vector *dataset, int *dataset_dim, int *dims, unsigned int min_fft, unsigned int max_fft, unsigned short num_iter);