#include<stdio.h>
#include<stdbool.h>

typedef struct level level;

struct level {
  unsigned int id;
  unsigned int num_cells;
  float cell_edge_length;
  bool is_leaf;
  cell * cells;
  struct level * prev; // previous level
  struct level * next; // next level
};

response init_first_level(pexeso_index * index);
response init_levels(pexeso_index * index);