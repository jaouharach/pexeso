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

void exit_with_error(char * message);
