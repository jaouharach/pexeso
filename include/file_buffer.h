#include<stdio.h>
#include <stdbool.h>
#include <stdlib.h>

struct file_buffer {
  struct cell * cell; // the buffer points back to its cell

  unsigned int disk_count; // number of vector in disk
  struct file_map * position_in_map; //the buffer points back to its position in file map
  
  v_type ** buffered_list; // list of vectors
  unsigned int buffered_list_size;   // number of vectors currently stored in this buffer

  bool in_disk; // false by default
  bool do_not_flush;
  
};

enum response file_buffer_init(struct cell * );
enum response flush_buffer_to_disk(struct pexeso_index *index, struct cell *cell);
enum response clear_file_buffer(struct pexeso_index *index, struct cell * cell);
enum response add_file_buffer_to_map(struct pexeso_index *index, struct cell *cell);
