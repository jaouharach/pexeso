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
typedef struct pairs pairs;
typedef struct candidate_pair candidate_pair;
typedef struct matching_pair matching_pair;
typedef struct query_result query_result;
typedef struct query_settings query_settings;
typedef struct inv_index inv_index;
typedef struct entry entry;
typedef struct best_fft best_fft;

// vector id (used to track vector ids of all vectors in one cell)
struct vid {
  unsigned int table_id;
  unsigned int set_id;
  unsigned int pos; // vector pos in set
  unsigned int set_size; // total vectors in set
  unsigned long set_pos_in_inv_index;
};

// set id 
struct sid {
  unsigned int table_id;
  unsigned int set_id;
  unsigned int set_size;
};

// vector
struct vector{
	unsigned int table_id;
	unsigned int set_id;
  unsigned int pos; // vector position in set
  unsigned int set_size;
  unsigned long set_pos_in_inv_index;
	v_type * values;
};

struct vector_tuple{
  struct vector * mtr_vector;
  struct vector * ps_vector;
};

// print error and exit program
void exit_with_failure(char * message);

// print warning
int warning(char *message);

// euclidean distance
float euclidean_distance(vector * v1, vector * v2, unsigned int v_len);

// euclidean distance for distance comparaison (without sqrt)
float euclidean_distance_cmp(vector * v1, vector * v2, unsigned int v_len);

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

/* get number of digits in integer */
int get_ndigits(unsigned int n);

/* print a progress bar */
void print_progress(int i, int j);

// compare set id returns 1 if set_id2 is greater than set_id1
int set_id_cmp(struct sid * set_id1, struct sid * set_id2);

// delete non empty directory
int remove_directory(const char *path);