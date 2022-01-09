#include <stdio.h>
#include <stdbool.h>

struct file_buffer_manager
{
  struct file_map *file_map;
  struct file_map *file_map_tail;
  int file_map_size;

  // memory array for metric space vectors
  char *mtr_memory_array;
  char *current_mtr_record; // next vector to be loaded in mtr_memory_array
  int current_mtr_record_index; 
  int max_mtr_record_index;

  // memory array for pivot space vectors
  char *ps_memory_array;
  char *current_ps_record; // next vector to be loaded in ps_memory_array
  int current_ps_record_index; 
  int max_ps_record_index;
};

struct file_map
{
  struct file_buffer *file_buffer;
  struct file_map *next;
  struct file_map *prev;
};

enum response init_file_buffer_manager(struct grid *grid);
enum response set_buffered_memory_size(struct grid * grid);
response get_file_buffer(struct grid *grid, struct cell *cell);