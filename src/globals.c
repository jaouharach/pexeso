#include <stdio.h>
#include <errno.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <ftw.h>
#include <sys/stat.h>
#include <unistd.h>
#include "../include/globals.h"

float euclidean_distance(vector * v1, vector * v2, unsigned int v_len)
{
    if (v1 == NULL || v2 == NULL)
        exit_with_failure("Error in globals.c: NULL pointer to vector.");

    if(isnan(v1->values[0]) || isnan(v1->values[0]))
        exit_with_failure("Error in globals.c: Cannot compute euclidean distance for nan vector.");

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

/* print warning and ask for user action, return 0 if No or 1 if user wants to continue */
int warning(char *message)
{
    fprintf(stderr, "%s (y/n):\n", message);
    char * resp = "n"; 
    scanf("%c", resp);
    if (resp[0] == 'n' || resp[0] == 'N')
        return 0;
}

void print_vector(vector * v, unsigned int v_len)
{
    printf("\t\t(");
    for(int i = 0; i < v_len; i++)
    {
        if(i == v_len - 1)
        {
            printf("%f", v->values[i]);
            break;
        }
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
    int * bin_arr;

    // allocate memory for result
    vector * result = malloc(sizeof(struct vector)* n);
    for (int i = 0; i < n; i++)
        result[i].values = malloc(sizeof(int) * dim);
    
    // cartesian product arr x arr of dimention k
    for (int i = 0; i < n; i++)
    {
        bin_arr = integer_to_binary_array(i, dim);
        for(int j = 0; j < dim; j++)
        {
            result[i].values[j] = (float) array[bin_arr[j]];
        }
        // free memory
        free(bin_arr);
    }

    return result;
}


/* add integer to array if it doesn't exist */
bool array_add(int * arr, int curr_size, int number)
{
    int i;
    for(i = 0; i < curr_size; i++)
    {
        if(number == arr[i])
            return false;
        
    }
    arr[i] = number; //add element to end of array
    return true;
}

/* get current time */
void get_current_time(char * time_buf)
{
    time_t timer;
    
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(time_buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);
}

/* get mean of a vector */
v_type get_vector_mean(vector * vector, unsigned int v_len)
{
    v_type mean = 0.0;
    for(int i = 0; i < v_len; i++)
    {
        mean += vector->values[i];
    }
    mean /= v_len;

    return mean;
}

/* get magnitude of a vector */
v_type get_vector_magnitude(vector * vector, unsigned int v_len)
{
    v_type mag = 0.0;

    for(int i = 0; i < v_len; i++)
    {
        mag += pow(vector->values[i], 2);
    }
    mag = sqrt(mag);

    return mag;
}

/* create directory, if directory exists ask user for action */
enum response create_grid_dir(const char * dir_path)
{
    struct stat sb;
    if (stat(dir_path, &sb) == 0 && S_ISDIR(sb.st_mode)) // if directory exists ask user for action
    {
        fprintf(stderr, "Grid root directory '%s' already exists!\n", dir_path);
        return FAILED;
    }

    mkdir(dir_path, 0777);
    return OK;
}

/* compare two pointers, returns true if both pointers point to the same object */
bool pointer_cmp(void * p1, void * p2)
{
    if
    (
        (p1 <= p2 && p1 >= p2) == true
        &&
        (p1 < p2) == false 
        &&
        (p1 > p2) == false
    )
    return true;

    return false;
}

/* get number of digits in integer */
int get_ndigits(unsigned int n)
{
	int total_digits = 0;
    if (n < 0)
        total_digits++;
	while (n != 0)
    {
		//4
		n = n/10;
		++total_digits;
	}
	return total_digits;
}