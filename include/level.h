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

response init_leaf_level(index_settings *, level *);
