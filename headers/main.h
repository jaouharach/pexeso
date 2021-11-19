#include<stdio.h>

char dataset_directory [] = "/home/jaouhara/Projects/pexeso/data_vectors.bin";

#define OK 1;
#define FAILED 0;

typedef struct level level;
typedef struct index index;
typedef struct index_settings index_settings;
typedef struct vector vector;
typedef struct cell cell;

struct index_settings {
    const char * root_directory;  
    unsigned int num_dim;
    float max_coordinate;    
    float min_coordinate;
    float leaf_cell_length;
};

struct index{
  unsigned long long total_records;
  level * first_level;
  struct index_settings * settings;
};

struct vector{
	int table_id;
	int set_id;
	float * values;
};

struct cell {
    struct cell * parent;
    char * filepath; //l[level]_center([center vector])_len[cell length]_size[num vectors].bin
    unsigned int is_leaf;
    vector * center;
    float length;
    unsigned long num_vectors;
};

struct level {
  unsigned int id;
  unsigned int num_cells;
  float cell_length;
  cell * cells;
  struct level * next_level;
};

//In progress
void insert_vector(char * index_directory, vector *);
void create_index(char * index_directory, char * dataset_directory);

void append_vector_to_cell(cell * , vector *);

// done.
int init_leaf_cells(level * leaf_level, index_settings * settings);
int init_leaf_cell(cell *, float length);

int init_leaf_level(index_settings *, level *);
int init_index(const char * root_directory,
                unsigned int num_dim,
                float max_coordinate,    
                float min_coordinate,
                float leaf_cell_length,
                index * pexeso_index);

void exit_with_error(char * message);


// file for each cell