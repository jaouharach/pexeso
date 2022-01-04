#include <stdio.h>
#include "../include/index.h"
#include "../include/file_buffer_manager.h"
#include "../include/file_buffer.h"

enum response init_file_buffer_manager(struct pexeso_index *index)
{
    index->buffer_manager = NULL;

    index->buffer_manager = malloc(sizeof(struct file_buffer_manager));
    if (index->buffer_manager == NULL)
        exit_with_failure("Error in file_buffer_manager.c: Could not allocate memory for buffer manager.\n");

    index->buffer_manager->max_buffered_size = 1000 * 1000 * 100;
    index->buffer_manager->current_vector_count = 0;

    index->buffer_manager->file_map = NULL;      //the file map is empty initially
    index->buffer_manager->file_map_tail = NULL; //the file map is empty initially

    index->buffer_manager->file_map_size = 0; // initially, zero files in file map

    // allocate memory for all vectors that will be indexed
    set_buffered_memory_size(index);

    return OK;
}

enum response set_buffered_memory_size(struct pexeso_index *index)
{
    if (index == NULL)
        exit_with_failure("Error in file_buffer_manager.c: Cannot set the buffered memory for a NULL index.");
    
    unsigned long buffered_memory_size_in_bytes = index->settings->buffered_memory_size * 1024 * 1024;
    index->buffer_manager->max_buffered_size = (long) (buffered_memory_size_in_bytes / sizeof(v_type));
    
    unsigned int vector_size_in_bytes = sizeof(v_type) * index->settings->mtr_vector_length;
    
    // why multiply by 2 !!! (I didn't here)
    unsigned long num_leaf_buffers =  buffered_memory_size_in_bytes /
                                    (unsigned long) (index->settings->max_leaf_size * vector_size_in_bytes);
    unsigned long leaf_buffer_size = sizeof(struct file_buffer) + (sizeof(v_type*) * index->settings->max_leaf_size);
    
    // total vector that can be stored in memory array = mem_array_size
    long long mem_array_size =  (long long)((buffered_memory_size_in_bytes -
                                            leaf_buffer_size * num_leaf_buffers) /
                                            vector_size_in_bytes);

    index->buffer_manager->memory_array = calloc(mem_array_size, vector_size_in_bytes);
        if(index->buffer_manager->memory_array == NULL)
            exit_with_failure("Error in file_buffer_manager.c:"
  		                        "Cannot allocate the requested buffer size.\n");

    index->buffer_manager->current_record_index = 0;
	index->buffer_manager->max_record_index = mem_array_size;
	index->buffer_manager->current_record = index->buffer_manager->memory_array;
	
    // printf("buffered_memory_size_in_bytes  = %lu\n", buffered_memory_size_in_bytes);
    // printf("max_buffered_size  = %llu\n", index->buffer_manager->max_buffered_size);
    // printf("vector_size_in_bytes  = %u\n", vector_size_in_bytes);
    // printf("num_leaf_buffers  = %lu\n", num_leaf_buffers);
    // printf("num_leaf_buffers  = %lu, leaf_buffer_size = %lu\n", num_leaf_buffers, leaf_buffer_size);
    // printf("mem_array_size  = %lld vectors\n", mem_array_size);
    // exit(1);
    
	return OK;
}