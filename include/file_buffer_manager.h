#include <stdio.h>
#include <stdbool.h>

struct file_buffer_manager
{
  struct file_map *file_map;
  struct file_map *file_map_tail;
  int file_map_size;

  unsigned long long max_buffered_size; // max number of v_type values that can be hold in memory in all file buffers
  long long current_values_count; // number of vectors currently in memory
  //   long batch_remove_size;

  char *memory_array;
  char *current_record; // next vector to be loaded in memory_array
  int current_record_index; 
  int max_record_index;
};

struct file_map
{
  struct file_buffer *file_buffer;
  struct file_map *next;
  struct file_map *prev;
};

enum response init_file_buffer_manager(struct pexeso_index *index);
enum response set_buffered_memory_size(struct pexeso_index * index);
response get_file_buffer(struct pexeso_index *index, struct cell *cell);