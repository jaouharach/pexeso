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
    unsigned int vector_size;
    float max_dim_value;    
    float min_dim_value;
    float leaf_cell_area;
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
  float cell_area;
  cell * cells;
  struct level * next_level;
};

//In progress
void insert_vector(char * index_directory, vector *v);
void create_index(char * index_directory, char * dataset_directory);
int init_cell(index_settings * settings, cell * c);
int init_level(index_settings * settings, level * l);
void append_vector_to_cell(cell * c, vector * v);

// done.
int init_index(const char * root_directory,
                unsigned int vector_size,
                float max_dim_value,    
                float min_dim_value,
                float leaf_cell_area,
                index * pexeso_index);

void exit_with_error(char * message);