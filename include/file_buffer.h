#include<stdio.h>
#include <stdbool.h>


struct file_buffer {

  cell * cell; // the buffer points back to its cell
  vector * buffered_list; // list of vectors

  int buffered_list_size;   // number of vectors currently stored in this buffer

  bool in_disk; // false by default
  bool do_not_flush;
  
};

response file_buffer_init(cell * );
response flush_buffer_to_disk(file_buffer *);
