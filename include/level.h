#include<stdio.h>

typedef struct level level;

struct level {
  unsigned int id;
  unsigned int num_cells;
  float cell_length;
  cell * cells;
  struct level * next_level;
};

response init_leaf_level(index_settings *, level *);
