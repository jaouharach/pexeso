#include <stdio.h>
#include "../include/index.h"
#include "../include/cell.h"
#include "../include/file_buffer.h"

response file_buffer_init(cell *cell)
{
    cell->file_buffer = NULL;
    cell->file_buffer = malloc(sizeof(struct file_buffer));

    if (cell->file_buffer == NULL)
        exit_with_error("Error in file_buffer.c: Could not allocate memory for file buffer.\n");

    cell->file_buffer->in_disk = false;
    
    cell->file_buffer->buffered_list = NULL;
    cell->file_buffer->buffered_list_size = 0;

    cell->file_buffer->cell = cell;
    cell->file_buffer->do_not_flush = false;

    return OK;
}


