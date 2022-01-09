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

typedef struct grid grid;
typedef struct grid_settings grid_settings;
typedef struct level level;
typedef struct cell cell;
typedef struct vector vector;
typedef struct vector_tuple vector_tuple;
typedef struct file_buffer_manager file_buffer_manager;
typedef struct file_buffer file_buffer;
typedef struct file_map file_map;
typedef struct candidate_pair candidate_pair;
typedef struct matching_pair matching_pair;
typedef struct query_result query_result;
typedef struct query_settings query_settings;
typedef struct inv_index inv_index;
typedef struct entry entry;

// vector id (used to track vector ids of all vectors in one cell)
struct vid {
  unsigned int table_id;
  unsigned int set_pos;
};
// vid (vector id) and sid (set id) are the same 
struct sid {
  unsigned int table_id;
  unsigned int set_pos;
};

// vector
struct vector{
	unsigned int table_id;
	unsigned int set_id;
	v_type * values;
};

struct vector_tuple{
	vector * mtr_vector;
  vector * ps_vector;
};

// print error and exit program
void exit_with_failure(char * message);

// print warning
int warning(char *message);

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

/* get mean of a vector */
v_type get_vector_mean(vector * vector, unsigned int v_len);

/* get magnitude of a vector */
v_type get_vector_magnitude(vector * vector, unsigned int v_len);

/* create directory, if directory exists ask user for action */
enum response create_grid_dir(const char * dir_path);

/* compare two pointers, returns true if both pointers point to the same object */
bool pointer_cmp(void * p1, void * p2);