#include <stdio.h>
#include <errno.h>
#include <float.h>
#include <stdlib.h>

#include "../include/globals.h"

float euclidean_distance(vector * v1, vector * v2, unsigned int v_len)
{
    if (v1 == NULL || v2 == NULL)
        exit_with_failure("Error in globals.c: NULL pointer to vector.");

    float d = 0.0;
    for(int j = 0; j < v_len; j++)
    {
        d = d + (v1->values[j] - v2->values[j]) * (v1->values[j] - v2->values[j]);
    }
    return (float) sqrt(d);
}

/* tranform vector to pivot space */
void transform_vector(vector * v, unsigned int v_len, vector * v_transf, vector * pivots, unsigned int num_pivots)
{
    for(int i = 0; i < num_pivots; i++)
        v_transf->values[i] = euclidean_distance(v, &pivots[i], v_len);
}

/* print error and kill process */
void exit_with_failure(char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

/* print warning and ask for user action */
void warning(char *message)
{
    fprintf(stderr, "%s\n", message);
    printf("Would you like to continue (y/n)? : ");
    char * resp = "n"; 
    scanf("%c", resp);
    if (resp[0] == 'n' || resp[0] == 'N')
        exit(1);
}

void print_vector(vector * v, unsigned int v_len)
{
    printf("(");
    for(int i = 0; i < v_len; i++)
    {
        printf("%f, ", v->values[i]);
    }
    printf(")\n");
}

/* copy v2 values in v1 */
void vector_cpy(vector * dest, vector * src, unsigned int v_len)
{
    for(int i = 0; i < v_len; i++)
    {
        dest->values[i] = src->values[i];
    }
}

/* convert integer to decimal array */
int * integer_to_binary_array(int number, unsigned int arr_len)
{
    int * bit_array = calloc(arr_len, sizeof(int));
    for(int i = arr_len - 1 ; number > 0 ; i--)
    {
        bit_array[i] = number%2;    
        number = number/2; 
    }
    return bit_array;
}

/* cartesian product (array x array) of dimension dim (array must be of length 2) */
vector * self_cartesian_product(int * array, unsigned int dim)
{
    // number of results
    int n = (int) pow(2, dim);
    // allocate memory for result
    vector * result = malloc(sizeof(struct vector)* n);
    for (int i = 0; i < n; i++)
        result[i].values = malloc(sizeof(int) * dim);
    
    // cartesian product arr x arr of dimention k
    for (int i = 0; i < n; i++)
    {
        int * bin_arr = integer_to_binary_array(i, dim);
        for(int j = 0; j < dim; j++)
        {
            result[i].values[j] = (float) array[bin_arr[j]];
        }
    }

    return result;
}


