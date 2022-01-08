#include<stdio.h>
#include<stdbool.h>

typedef struct level level;

struct level {
  unsigned int id;
  unsigned int num_cells;
  float cell_edge_length;
  bool is_leaf;
  bool is_first;
  bool is_root;
  struct cell * cells;
  struct level * next; // next level
};

/* initialize root level */
enum response init_root(struct grid *index);
enum response init_first_level(struct grid * index);
enum response init_levels(struct grid * index);