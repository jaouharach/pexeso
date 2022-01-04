#include<stdio.h>
#include <stdbool.h>


struct file_buffer {
  cell * cell; // the buffer points back to its cell

  unsigned int disk_count;
  struct dstree_file_map * position_in_map; //the buffer points back to its position in file map
  
  vector ** buffered_list; // list of vectors
  unsigned int buffered_list_size;   // number of vectors currently stored in this buffer

  bool in_disk; // false by default
  bool do_not_flush;
  
};

response file_buffer_init(cell * );
response flush_buffer_to_disk(file_buffer *);
