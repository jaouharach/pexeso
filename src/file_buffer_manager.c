#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/hgrid.h"
#include "../include/cell.h"
#include "../include/file_buffer_manager.h"
#include "../include/file_buffer.h"

enum response init_file_buffer_manager(struct grid *grid)
{
    grid->buffer_manager = NULL;

    grid->buffer_manager = malloc(sizeof(struct file_buffer_manager));
    if (grid->buffer_manager == NULL)
        exit_with_failure("Error in file_buffer_manager.c: Could not allocate memory for buffer manager.\n");

    grid->buffer_manager->file_map = NULL;      //the file map is empty initially
    grid->buffer_manager->file_map_tail = NULL; //the file map is empty initially

    grid->buffer_manager->file_map_size = 0; // initially, zero files in file map

    // allocate memory for all vectors that will be indexed
    set_buffered_memory_size(grid);

    return OK;
}

enum response set_buffered_memory_size(struct grid *grid)
{
    if (grid == NULL)
        exit_with_failure("Error in file_buffer_manager.c: Cannot set the buffered memory for a NULL grid.");

    /*
        memory used by the program is muchlarger than buffered_memory_size
        in this version buffered_memory_size is used for string vectors only and it does not 
        take into account memory required to store structs (file_buffer, cell, level ...)
    */
    unsigned long num_mtr_bytes = grid->settings->mtr_buffered_memory_size * 1024 * 1024;
    unsigned long num_ps_bytes = grid->settings->ps_buffered_memory_size * 1024 * 1024;
    
    unsigned long mtr_vector_size_in_bytes = sizeof(v_type) * grid->settings->mtr_vector_length;
    unsigned long ps_vector_size_in_bytes = sizeof(v_type) * grid->settings->num_pivots;

    long max_mtr_buffered_vectors = (long)(num_mtr_bytes / mtr_vector_size_in_bytes);
    long max_ps_buffered_vectors = (long)(num_ps_bytes / ps_vector_size_in_bytes);

    printf("max_mtr_buffered_vectors  = %lu\n", max_mtr_buffered_vectors);
    printf("max_ps_buffered_vectors  = %lu\n", max_ps_buffered_vectors);


    grid->buffer_manager->mtr_memory_array = calloc(max_mtr_buffered_vectors, mtr_vector_size_in_bytes);
    if (grid->buffer_manager->mtr_memory_array == NULL)
        exit_with_failure("Error in file_buffer_manager.c:"
                          " Cannot allocate the requested buffer size.\n");

    grid->buffer_manager->ps_memory_array = calloc(max_ps_buffered_vectors, ps_vector_size_in_bytes);
    if (grid->buffer_manager->ps_memory_array == NULL)
        exit_with_failure("Error in file_buffer_manager.c:"
                          " Cannot allocate the requested buffer size.\n");

    grid->buffer_manager->current_ps_record_index = 0;
    grid->buffer_manager->max_mtr_record_index = max_mtr_buffered_vectors;
    grid->buffer_manager->current_mtr_record = grid->buffer_manager->mtr_memory_array;

    grid->buffer_manager->current_mtr_record_index = 0;
    grid->buffer_manager->max_ps_record_index = max_ps_buffered_vectors;
    grid->buffer_manager->current_ps_record = grid->buffer_manager->ps_memory_array;

    return OK;
}

response get_file_buffer(struct grid *grid, struct cell *cell)
{
    if (cell->file_buffer == NULL) // cell does not have a file buffer yet
    {
        if (!file_buffer_init(cell))
            exit_with_failure("Error in file_buffer_manager.c: Could not initialize the file buffer \
                           for the node.\n");

        if(cell->file_buffer == NULL)
            exit_with_failure("Error in file_buffer_manager.c: Could not initialize the file buffer \
                           for the node.\n");

        if (!add_file_buffer_to_map(grid, cell))
            exit_with_failure("Error in file_buffer_manager.c: Could not add the file buffer \
                           to the map.\n");

        if(cell->file_buffer == NULL)
            exit_with_failure("Error in file_buffer_manager.c: Could not initialize the file buffer \
                           for the node.\n");
    }

    // if buffer limit has been reached (todo) flush huffer to disk
    int mtr_buffer_limit = grid->buffer_manager->max_mtr_record_index;
    int ps_buffer_limit = grid->buffer_manager->max_ps_record_index;

    if (grid->buffer_manager->current_mtr_record_index > mtr_buffer_limit)
        exit_with_failure("Error in file_buffer_manager.c: memory array has been completly filled! cannot add any more vectors. please increase  buffer memory size.");
    
    if (grid->buffer_manager->current_ps_record_index > ps_buffer_limit)
        exit_with_failure("Error in file_buffer_manager.c: memory array has been completly filled! cannot add any more vectors. please increase  buffer memory size.");
    
    // if (grid->buffer_manager->current_record_index > buffer_limit)
    // {
    //     char *curr_time;
    //     curr_time = NULL;
    //     curr_time = malloc(sizeof(char) * 26);
    //     get_current_time(curr_time);

    //     printf("%s, batch remove ! %u  \n", curr_time, grid->buffer_manager->current_record_index);
    //     free(curr_time);

    //     unsigned long to_size = 0;
    //     struct file_map *currP = NULL;
    //     currP = grid->buffer_manager->file_map;

    //     while (currP != NULL)
    //     {
    //         //flush buffer with position idx in the file map of this grid
    //         if (!flush_buffer_to_disk(grid, currP->file_buffer->cell)) //flush the actual buffer of the node
    //             exit_with_failure("Error in file_buffer_manager.c: Could not flush the buffer \
    //                                 for node to disk.\n");

    //         currP = currP->next;
    //     }
    //     memset(grid->buffer_manager->memory_array, 0, grid->buffer_manager->max_record_index);
    //     grid->buffer_manager->current_record_index = 0;
    //     grid->buffer_manager->current_record = grid->buffer_manager->memory_array;
    // }

    return OK;
}