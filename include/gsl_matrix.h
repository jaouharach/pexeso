#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

/* inverse of a square gsl_matrix */
gsl_matrix * gsl_matrix_inverse(gsl_matrix * input_matrix, int *dim);

/* matrix x matrix multiplication */
gsl_matrix * gsl_matrix_product(gsl_matrix *matrix_a, int *dim_a, gsl_matrix *matrix_b, int *dim_b);

/* transpose of a data_set matrix */
gsl_matrix * gsl_matrix_get_transpose(gsl_matrix * input_matrix, int *dim);
void gsl_matrix_transpose_in_place(gsl_matrix * input_matrix, int *dim);

/* get determinant of a square matrix (dim[0] number of rows = dim[1] number of columns) */ 
double gsl_matrix_determinant(gsl_matrix *matrix_a, int dim);

/* get rank of a matrix, rank = number of linearly independant rows */ 
unsigned int gsl_matrix_rank(gsl_matrix *matrix, int *dim);

/* get covariance matrix of a matrix */
gsl_matrix * gsl_matrix_covariance(gsl_matrix * input_matrix, int *dim);

/* Mean centerize a gsl matrix */
void gsl_matrix_centerize(gsl_matrix *input_matrix, int *dim);

/* orthonormalization of a matrix */
gsl_matrix * orthonormalization(gsl_matrix *matrix, int *dim, int is_large_data_set);

/* change values in a (h x w) sub matrix of the matrix starting at [row][col] */
void gsl_matrix_set_part(gsl_matrix *matrix, unsigned int row, unsigned int col, gsl_matrix * changes, unsigned int height, unsigned int width);

/* get a (h x w) sub matrix of the matrix starting at [row][col] */
gsl_matrix * gsl_matrix_get_part(gsl_matrix *matrix, unsigned int row, unsigned int col, unsigned int height, unsigned int width);

/* flip columns of a matrix */
gsl_matrix * gsl_matrix_flip_columns(gsl_matrix *matrix, int *dim);
/* flip vector */
gsl_vector * vector_flip(gsl_vector * v, int dim);
/* swap rows in matrix */
void gsl_rows_swap(gsl_matrix *matrix, int row_1, int row_2, int num_cols);

/* convert matrix to positive matrix */
void gsl_to_positive_matrix(gsl_matrix *input_matrix, int *dim);

/* sort the matrix rows into ascending order, according to the natural ordering of the matrix values in the given column. */
gsl_matrix * gsl_matrix_sort_by_column(gsl_matrix * input_matrix, int *dim, int col);

/* print gsl_matrix */
void gsl_matrix_print(gsl_matrix * matrix, int * dim);

/* convert matrix from type (vector *) to type gsl_matrix */
gsl_matrix * to_gsl_matrix(vector * matrix, int *dim);
/* convert matrix from type gsl_matrix to type (vector *) */
vector * to_original_matrix(gsl_matrix * matrix, int *dim);