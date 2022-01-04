#include "../include/index.h"
#include "../include/gsl_matrix.h"


/* inverse of a square gsl_matrix */
gsl_matrix * gsl_matrix_inverse(gsl_matrix * input_matrix, int *dim)
{
    gsl_permutation *p = gsl_permutation_alloc(dim[0]);
    gsl_matrix *inverse_matrix = gsl_matrix_alloc(dim[0], dim[0]);
    int s;
    
    gsl_linalg_LU_decomp(input_matrix, p, &s); // Compute the LU decomposition of this matrix
    gsl_linalg_LU_invert(input_matrix, p, inverse_matrix); // Compute the  inverse of the LU decomposition

    gsl_permutation_free(p);

    return inverse_matrix;

}
/* matrix x matrix multiplication */
gsl_matrix * gsl_matrix_product(gsl_matrix *matrix_a, int *dim_a, gsl_matrix *matrix_b, int *dim_b)
{
    if (dim_a[1] != dim_b[0])
        exit_with_failure("Error in globals.c: Couldn't compute matrix product, number of columns in first matrix must be equal to number of rows in the second matrix!");

    // dimension of the product matrix = num_rows_in_a x num_cols_in_b
    int num_rows = dim_a[0];
    int num_cols = dim_b[1];
    // int dim_result [] = {num_rows, num_cols};

    gsl_matrix *result = gsl_matrix_alloc(num_rows, num_cols);
    if (result == NULL)
        exit_with_failure("Error in globals.c: Couldn't allocate memory for matrix product.");

    // matrix x matrix multiplication
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, matrix_a, matrix_b, 0.0, result);

    return result;
}
/* transpose of a data_set matrix */
gsl_matrix * gsl_matrix_get_transpose(gsl_matrix * input_matrix, int *dim)
{
    int num_rows = dim[0];
    int num_cols = dim[1];
    gsl_matrix *t = gsl_matrix_alloc(num_cols, num_rows);
    if(t == NULL)
        exit_with_failure("Error in globals.c: Couldn't allocate memory for transpose matrix.");

    for (int i = 0; i < num_rows; i++)
        for (int j = 0; j < num_cols; j++)
            gsl_matrix_set(t, j, i, gsl_matrix_get(input_matrix, i, j));
        
    return t;
}
void gsl_matrix_transpose_in_place(gsl_matrix * input_matrix, int *dim)
{
    int num_rows = dim[0];
    int num_cols = dim[1];
    gsl_matrix *t = gsl_matrix_alloc(num_cols, num_rows);
    if(t == NULL)
        exit_with_failure("Error in globals.c: Couldn't allocate memory for transpose matrix.");

    for (int i = 0; i < num_rows; i++)
        for (int j = 0; j < num_cols; j++)
            gsl_matrix_set(t, j, i, gsl_matrix_get(input_matrix, i, j));
    
    input_matrix = t;
    gsl_matrix_free(t);
}
/* get determinant of a square matrix (dim[0] number of rows = dim[1] number of columns) */ 
double gsl_matrix_determinant(gsl_matrix *matrix_a, int dim)
{
    gsl_matrix *matrix_b = gsl_matrix_alloc(dim, dim);
    if (matrix_b == NULL)
        exit_with_failure("Error in globals.c: Could't allocate memory for matrix in determinant() function.");
   
    float s = 1, det = 0;
    int i, j, m, n, c;

    if (dim == 1)
        return gsl_matrix_get(matrix_a, 0, 0);

    det = 0;
    for (c = 0; c < dim; c++)
    {
        m = 0;
        n = 0;
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
            {
                gsl_matrix_set(matrix_b, i, j, 0);
                if (i != 0 && j != c)
                {
                    gsl_matrix_set(matrix_b, m, n, gsl_matrix_get(matrix_a, i, j));
                    if (n < (dim - 2))
                        n++;
                    else
                    {
                        n = 0;
                        m++;
                    }
                }
            }
        }
        det = det + s * (gsl_matrix_get(matrix_a, 0, c) * gsl_matrix_determinant(matrix_b, dim - 1));
        s = -1 * s;
    }
    printf("det = %f\n", det);
    return det;
}
/* get rank of a matrix, rank = number of linearly independant rows */ 
unsigned int gsl_matrix_rank(gsl_matrix *matrix, int *dim)
{
    unsigned int num_rows = dim[0];
    unsigned int num_cols = dim[1];
    unsigned int rank = num_cols; // initially rank = num_cols 
    
    //make a copy of the input matrix
    gsl_matrix *matrix_copy = gsl_matrix_alloc(num_rows, num_cols);
    if(matrix_copy == NULL)
        exit_with_failure("Error int globals.c: Couldn't allocate memory for matrix copy in function gsl_matrix_rank().");
    gsl_matrix_memcpy(matrix_copy, matrix);

    for (int r = 0; r < rank; r++)
    {
        // check that diagonal element is not zero
		if (gsl_matrix_get(matrix_copy, r, r) != 0)
		{
            for (int c = 0; c < num_rows; c++)
            {
                if (c != r)
                {
                    // This makes all entries of current column as 0 except diognal entry matrix[r][r]
                    double mult = (double) gsl_matrix_get(matrix_copy, c, r) / gsl_matrix_get(matrix_copy, r, r);
                    for (int i = 0; i < rank; i++)
                    {
                        double val = gsl_matrix_get(matrix_copy, c, i) - (mult * gsl_matrix_get(matrix_copy, r, i));
                        gsl_matrix_set(matrix_copy, c, i,  val);
                    }
                        
                }
            }
		}
        /*  if diognal element is already zero. Two cases
            arise:
                1) If there is a row below it with non-zeroentry, then swap this row with that row and process that row
                2) If all elements in current column below mat[r][row] are 0, then remove this column by swapping it with 
                    last column and reducing number of columns by 1. 
        */
        else
        {
			bool reduce = true;

			/* Find the non-zero element in current column */
			for (int i = r + 1; i < num_rows; i++)
			{
				// Swap the row with non-zero element with this row.
				if (gsl_matrix_get(matrix_copy, i, r) != 0)
				{
					gsl_rows_swap(matrix_copy, r, i, rank);
					reduce = false;
					break;
				}
			}
			/*  
                If we did not find any row with non-zero element in current column, then all
			    values in this column are 0.
            */
			if (reduce)
			{
				// Reduce number of columns
				rank--;

				// Copy the last column here
				for (int i = 0; i < num_rows; i++)
                    gsl_matrix_set(matrix_copy, i, r,  gsl_matrix_get(matrix_copy, i, rank));
			}
			// Process this row again
			r--;
		}
    }

    gsl_matrix_free(matrix_copy);
    return rank;
}
/* get covariance matrix of a matrix */
gsl_matrix * gsl_matrix_covariance(gsl_matrix * input_matrix, int *dim)
{
    unsigned int num_cols = dim[1];

    /* Allocate memory for covariance matrix */
    gsl_matrix * matrix_cov = gsl_matrix_alloc(num_cols, num_cols);
    
    double cov;
    for (int i = 0; i < num_cols; i++) {
        _gsl_vector_view a = gsl_matrix_column(input_matrix, i);
        for (int j = 0; j < num_cols; j++) {
            _gsl_vector_view b = gsl_matrix_column(input_matrix, j);
            cov = gsl_stats_covariance(a.vector.data, a.vector.stride,b.vector.data, b.vector.stride, a.vector.size);
            gsl_matrix_set (matrix_cov, i, j, cov);
        }
    }
    
    return matrix_cov;
}

/* Mean centerize a gsl matrix */
void gsl_matrix_centerize(gsl_matrix *input_matrix, int *dim)
{
    int num_rows = dim[0];
    int num_cols = dim[1];

    gsl_vector  *means = gsl_vector_alloc(num_cols);
    if (means == NULL)
        exit_with_failure("Error in globals.c: Couldn't allocate memory for mean matrix.");

    for (int i = 0; i < num_cols; i++)
    {
        for (int j = 0; j < num_rows; j++)
            gsl_vector_set(means, i, gsl_matrix_get(input_matrix, j, i) + gsl_vector_get(means, i));

        gsl_vector_set(means, i, gsl_vector_get(means, i) / num_rows);
    }

    for (int i = 0; i < num_cols; i++)
    {
        for (int j = 0; j < num_rows; j++)
            gsl_matrix_set(input_matrix, j, i, gsl_matrix_get(input_matrix, j, i) - gsl_vector_get(means, i));
    }

    gsl_vector_free(means);
}

/* orthonormalization of a matrix */
gsl_matrix * orthonormalization(gsl_matrix *matrix_a, int *dim, int is_large_data_set)
{
    unsigned int num_rows = dim[0];
    unsigned int num_cols = dim[1];
    unsigned int rank;


    // printf("orthonormalization:\n C Matrix:\n");
    // gsl_matrix_print(matrix_a, dim);

    // sigular value decompostion of input matrix
    // printf("Run SVD...\n");
    gsl_matrix * U = matrix_a;
    gsl_vector * S = gsl_vector_alloc(num_cols);
    gsl_matrix * V = gsl_matrix_alloc(num_cols, num_cols);
    gsl_vector * work = gsl_vector_alloc(num_cols);
    gsl_matrix * X = gsl_matrix_alloc(num_cols, num_cols);


    if(!is_large_data_set) // num_rows ! >> num_cols
        gsl_linalg_SV_decomp(matrix_a, V, S, work);
    else // num_rows >> num_cols
        gsl_linalg_SV_decomp_mod(U, X, V, S, work);
    
    // rank of  matrix
    rank = gsl_matrix_rank(U, dim);
    // printf("Rank = %u\n", rank);

    gsl_matrix * result = gsl_matrix_alloc(num_rows, rank);
    gsl_vector * temp_column = gsl_vector_alloc(num_rows);
    gsl_vector * column_ri = gsl_vector_alloc(num_rows);

    // copy first column to result
    gsl_vector * u_col_0 = gsl_vector_alloc(num_rows);
    gsl_matrix_get_col(u_col_0, U, 0);
    gsl_matrix_set_col (result, 0, u_col_0);
    
    for(int c = 1; c < rank; c++)
    {
        gsl_matrix_get_col(column_ri, U, c);
        gsl_matrix_set_col(result, c, column_ri);

        for(int i = 0; i < c; i++)
        {
            /* compute: temp_column = matrix_u[][c] x result[][i] / norm2(result[][i]) */

            gsl_matrix_get_col(column_ri, result, i); // column_ri = resutl[][i]
            gsl_matrix_get_col(temp_column, U, i); // temp_column = matrix_u[][i]
            gsl_vector_mul(temp_column, column_ri); // temp_column = matrix_u[][c] x result[][i]
            double norm2 = gsl_blas_dnrm2(temp_column); // norm2(result[][i])
            gsl_vector_scale(temp_column, 1/norm2);


            /* for every column j in temp_column: temp_column[j] =  temp_column[j] *  result[j][i] */
            for(int j = 0; j < num_rows; j++)
                gsl_vector_set(temp_column, j, gsl_vector_get(temp_column, j) * gsl_matrix_get(result, j, i));

            /* for every value k in result[][c]: result[k][c] = result[k][c] - temp_column[k] */
            for(int k = 0; k < num_rows; k++)
                gsl_matrix_set(result, k, c, gsl_matrix_get(result, k, c) - gsl_vector_get(temp_column, k));
        }
    }

    for(int c = 0; c < rank; c++)
    {
        // temp_column = sqrt(norm2(result[][c]))
        gsl_matrix_get_col(temp_column, result, c); 
        double norm2 = gsl_blas_dnrm2(temp_column);
        double sqrt_norm2 = sqrt(norm2);
        gsl_vector_set_all(temp_column, sqrt_norm2);

        // for every row i in current column c: result[i][c] = result[i][c] / temp_column[i]
        for(int i = 0; i < num_rows; i++)
                gsl_matrix_set(result, i, c, gsl_matrix_get(result, i, c) / gsl_vector_get(temp_column, i));
    }


    // printf("U matrix:\n");
    // gsl_matrix_print(U, dim);

    // printf("Vector S:\n");
    // for(int i = 0; i < rank; i++)
    //     printf("\t%.3f", gsl_vector_get(S, i));

    // printf("matrix V:\n");
    // int dim_v [] = {num_cols, num_cols};
    // gsl_matrix_print(V, dim_v);

    gsl_vector_free(temp_column);
    gsl_vector_free(column_ri);
    gsl_vector_free(S);
    gsl_matrix_free(result);
    gsl_matrix_free(V);
    gsl_matrix_free(X);

    return U;
}
/* change values in a (h x w) sub matrix of the matrix starting at [row][col] */
void gsl_matrix_set_part(gsl_matrix *matrix, unsigned int row, unsigned int col, gsl_matrix * changes, unsigned int height, unsigned int width)
{
    for(int i = 0; i < height; i++)
    {
        for(int j = 0; j < width; j++)
            gsl_matrix_set(matrix, i+row, j+col, gsl_matrix_get(changes, i, j));

    }
}
/* get a (h x w) sub matrix of the matrix starting at [row][col] */
gsl_matrix * gsl_matrix_get_part(gsl_matrix *matrix, unsigned int row, unsigned int col, unsigned int height, unsigned int width)
{
    gsl_matrix * sub_matrix = gsl_matrix_alloc(height, width);
    if(sub_matrix == NULL)
        exit_with_failure("Error in gsl_matrix.c: Couldn't allocate memory for sub matrix.");

    for(int i = 0; i < height; i++)
    {
        for(int j = 0; j < width; j++)
            gsl_matrix_set(sub_matrix, i, j, gsl_matrix_get(matrix, i+row, j+col));

    }
    return sub_matrix;
}
/* flip columns of a matrix, example at: https://dst.lbl.gov/ACSSoftware/colt/api/cern/colt/matrix/DoubleMatrix2D.html#viewColumnFlip() */
gsl_matrix * gsl_matrix_flip_columns(gsl_matrix *matrix, int *dim)
{
    unsigned int num_rows = dim[0];
    unsigned int num_cols = dim[1];

    gsl_matrix * fliped_matrix = gsl_matrix_alloc(num_rows, num_cols);
    if(fliped_matrix == NULL)
        exit_with_failure("Error in globals.c: Couldn't allocate memory for fliped matrix in function matrix_flip_columns().");
    

    for(int i = 0; i < num_rows; i++)
    {
        for(int j = 0; j < num_cols; j++)
        {
            int k = num_cols - j - 1;
            gsl_matrix_set(fliped_matrix, i, j, gsl_matrix_get(matrix, i, k));
        }
    }

    return fliped_matrix;
}
/* flip vector */
gsl_vector * vector_flip(gsl_vector * v, int dim)
{
    gsl_vector * fliped_vector = gsl_vector_alloc(dim);
    if(fliped_vector == NULL)
        exit_with_failure("Error in globals.c: Couldn't allocate memory for fliped vector in function vector_flip().");

    for(int j = 0; j < dim; j++)
        gsl_vector_set(fliped_vector, j, gsl_vector_get(v, dim - j - 1));
    
    return fliped_vector;
}
/* swap two rows of a matrix */ 
void gsl_rows_swap(gsl_matrix *matrix, int row_1, int row_2, int num_cols)
{
    for (int i = 0; i < num_cols; i++)
    {
        double temp = gsl_matrix_get(matrix, row_1, i);
        gsl_matrix_set(matrix, row_1, i,  gsl_matrix_get(matrix, row_2, i));
        gsl_matrix_set(matrix, row_2, i,  temp);
    }
}

/* convert matrix to positive matrix */
void gsl_to_positive_matrix(gsl_matrix *input_matrix, int *dim)
{
    unsigned int num_rows = dim[0];
    unsigned int num_cols = dim[1];

    for(int i = 0; i < num_rows; i++)
        for(int j = 0; j < num_cols; j++)
            gsl_matrix_set(input_matrix, i, j, fabs(gsl_matrix_get(input_matrix, i, j)));
}

/* sort the matrix rows into ascending order, according to the natural ordering of the matrix values in the given column. */
gsl_matrix * gsl_matrix_sort_by_column(gsl_matrix * input_matrix, int *dim, int col)
{
    unsigned int num_rows = dim[0];
    unsigned int num_cols = dim[1];

    gsl_matrix * sorted_matrix = gsl_matrix_alloc(num_rows, num_cols);
    gsl_matrix_memcpy(sorted_matrix, input_matrix);

    double temp;

    for(int i = 0; i < num_rows; i++)
    {
            for(int j = i + 1;j < num_rows;j++)
            {
                    if(gsl_matrix_get(sorted_matrix, i, col) > gsl_matrix_get(sorted_matrix, j, col))
                    {
                        for(int x = 0; x < num_cols; x++) {
                                temp = gsl_matrix_get(sorted_matrix, i, x);
                                gsl_matrix_set(sorted_matrix, i, x,  gsl_matrix_get(sorted_matrix, j, x));
                                gsl_matrix_set(sorted_matrix, j, x,  temp);
                            }
                    }
            }
    }
    return sorted_matrix;
}

/* print gsl_matrix */
void gsl_matrix_print(gsl_matrix * matrix, int * dim)
{
    unsigned int num_rows = dim[0];
    unsigned int num_cols = dim[1];

    for(int i = 0; i < num_rows; i++)
    {
        for(int j = 0; j < num_cols; j++)
            printf("\t%.4f", gsl_matrix_get(matrix, i, j));
        printf("\n");
    }
}


/* convert matrix from type (vector *) to type gsl_matrix */
gsl_matrix * to_gsl_matrix(vector * input_matrix, int *dim)
{
    unsigned int num_rows = dim[0];
    unsigned int num_cols = dim[1];

    gsl_matrix * output_matrix = gsl_matrix_alloc(num_rows, num_cols);
    for(int i = 0; i < num_rows; i++)
        for(int j = 0; j < num_cols; j++)
            gsl_matrix_set(output_matrix, i, j, (double) input_matrix[i].values[j]);

    return output_matrix;
}
/* convert matrix from type gsl_matrix to type (vector *) */
vector * to_original_matrix(gsl_matrix * input_matrix, int *dim)
{
    unsigned int num_rows = dim[0];
    unsigned int num_cols = dim[1];

    vector * output_matrix = malloc(sizeof(struct vector) * num_rows);
    for(int i = 0; i < num_rows; i++)
    {
        output_matrix[i].values  = malloc(sizeof(v_type) * num_cols);
        for(int j = 0; j < num_cols; j++)
            output_matrix[i].values[j] = (v_type) gsl_matrix_get(input_matrix, i, j);
    }
    return output_matrix;
}
