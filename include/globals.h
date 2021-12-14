#include<stdio.h>

/* enumerations */
typedef enum response {OK = 1, FAILED = 0} response;


/* types */
typedef float v_type; // vector values

typedef struct pexeso_index pexeso_index;
typedef struct index_settings index_settings;
typedef struct level level;
typedef struct cell cell;
typedef struct vector vector;
typedef struct file_buffer file_buffer;

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

// FFT: farthest first traversal, n = data_set size (total vectors), k  = number of outliers to be retrieved
vector * fft(vector * data_set, unsigned int n, unsigned int k, unsigned int v_len);
float min_distance(vector * vectors, unsigned int num_vectors, vector * v, unsigned int v_len);
void print_vector(vector * v, unsigned int v_len);
void vector_cpy(vector * dest, vector * src, unsigned int v_len);