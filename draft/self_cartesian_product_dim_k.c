#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Find cartesian product of arr x arr of dimension k
typedef struct vector{
    int * values;
}vector;

void print_vector(vector vec, int dim);
int * to_binary_array(int input_int, int k);
int main()
{
    int arr[] = {1, -1};
    int r = 2;
    int k = 4;
    int n = (int) pow(2, k); // total number of vectors = 2^k
    
    // allocate memory for result
    vector * result = malloc(sizeof(struct vector)* n);
    for (int i = 0; i < n; i++)
        result[i].values = malloc(sizeof(int) * k);
    
    // cartesian product arr x arr of dimention k
    for (int i = 0; i < n; i++)
    {
        int * bin_arr = to_binary_array(i, k);
        for(int j = 0; j < k; j++)
        {
            result[i].values[j] = arr[bin_arr[j]];
        }
    }
    // print result
    for(int i = 0; i < n; i++)
        print_vector(result[i], k);
}

void print_vector(vector vec, int dim)
{
    printf("(");
    for(int j = 0; j < dim; j++)
    {
        printf("%d, ", vec.values[j]);
    }
    printf(")\n");
}

int * to_binary_array(int input_int, int k)
{
    int * bytearray = calloc(k, sizeof(int));
    for(int i = k-1 ; input_int > 0 ; i--)
    {
        bytearray[i] = input_int%2;    
        input_int = input_int/2; 
    }
    return bytearray;
}


// code equivalent in python:
/*
from itertools import product
arr = [1, -1]
k = 3
print (list(product(set(arr),repeat = k)))
*/
