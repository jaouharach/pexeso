#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <float.h>
#include <stdlib.h>
#include <stdbool.h>

/* enumerations */
typedef enum response {OK = 1, FAILED = 0} response;

/* types */
typedef float v_type; // vector values

typedef struct pexeso_index pexeso_index;
typedef struct index_settings index_settings;
typedef struct level level;
typedef struct cell cell;
typedef struct vector vector;
typedef struct file_buffer_manager file_buffer_manager;
typedef struct file_buffer file_buffer;
typedef struct file_map file_map;

struct vector{
	int table_id;
	int set_id;
	v_type * values;
};

// print error and exit program
void exit_with_failure(char * message);

// print warning
void warning(char *message);

// euclidean distance
float euclidean_distance(vector * v1, vector * v2, unsigned int v_len);

/* copy vector */
void vector_cpy(vector * dest, vector * src, unsigned int v_len);

/* cartesian product (array x array) of dimension dim  */
vector * self_cartesian_product(int *array, unsigned int dim);

/* convert integer to decimal array */
int * integer_to_binary_array(int number, unsigned int arr_len);

/* print matrix */
void matrix_print(vector * matrix, int * dim);

/* print vector */
void print_vector(vector * v, unsigned int v_len);

/* add integer to array if it doesn't exist */
bool array_add(int * arr, int curr_size, int number);

/* get current time */
void get_current_time(char * time_buf);