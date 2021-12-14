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
        d += (v1->values[j] - v2->values[j]) * (v1->values[j] - v2->values[j]);
    }
    return d;
}

// Farthest First traversal
float min_distance(vector * vectors, unsigned int num_vectors, vector * v, unsigned int v_len)
{
    float min_d = FLT_MAX, d;
    for (int i = 0; i < num_vectors; i++)
    {
        d = euclidean_distance(v, &vectors[i], v_len);
        if(d < min_d)
            min_d = d;
    }
    return min_d;
}

vector * fft(vector * data_set, unsigned int n, unsigned int k, unsigned int v_len)
{
    vector * outliers = malloc(sizeof(struct vector) * k);
    for(int i = 0; i < k; i++)
    {
        outliers[i].values = malloc(sizeof(v_type)* v_len);
    }
    // unsigned int  * outliers_idx = malloc(sizeof(unsigned int) * k);
    float bsf_d;
    int bsf_v;

    if (k == n)
        warning("Warning in globals.c: number of outliers is equal to total vectors in the dataset.");

    // pick an arbitrary vector as first outlier
    unsigned int rand_idx = (rand() % (n - 1));
    vector_cpy(&outliers[0], &data_set[rand_idx], v_len);


    // find k outliers, skip first
    for(int i = 1; i < k; i++)
    {
        bsf_d = FLT_MIN; // max distance to the current set of outliers
        bsf_v = -1; // position of the farthest vector in the data set.
        // find the farthest vector to the vectors in outliers list
        for(int j = 0; j < n; j++)
        {
            // check if vector is already in oultiers
            // (...)
            float d = min_distance(outliers, i+1, &data_set[j], v_len);
            if( d > bsf_d)
            {
                bsf_d = d;
                bsf_v = j;
            }
        }
        outliers[i].values = data_set[bsf_v].values;
    }
    return outliers;
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
        printf("%2f, ", v->values[i]);
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
