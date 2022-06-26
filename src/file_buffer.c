#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/hgrid.h"
#include "../include/cell.h"
#include "../include/file_buffer.h"
#include "../include/file_buffer_manager.h"
#include "../include/stats.h"

enum response file_buffer_init(struct cell *cell)
{
    cell->file_buffer = NULL;
    cell->file_buffer = (struct file_buffer *)malloc(sizeof(struct file_buffer));
    if (cell->file_buffer == NULL)
        exit_with_failure("Error in file_buffer.c: Could not allocate memory for file buffer.\n");

    cell->file_buffer->in_disk = false;
    cell->file_buffer->disk_count = 0; // total vectors in disk

    cell->file_buffer->mtr_buffered_list = NULL;
    cell->file_buffer->ps_buffered_list = NULL;

    cell->file_buffer->buffered_list_size = 0;

    cell->file_buffer->cell = cell;
    cell->file_buffer->position_in_map = NULL;

    cell->file_buffer->do_not_flush = false;

    return OK;
}

/* add file buffer to file map */
enum response add_file_buffer_to_map(struct grid *grid, struct cell *cell)
{
    // printf("Adding file_buffer of cell to map.\n");
    int idx = grid->buffer_manager->file_map_size;

    // first file
    if (idx == 0)
    {
        grid->buffer_manager->file_map = malloc(sizeof(struct file_map));
        if (grid->buffer_manager->file_map == NULL)
            exit_with_failure("Error in file_buffer_manager.c:"
                              "Could not allocate memory for the file map.\n");

        grid->buffer_manager->file_map[idx].file_buffer = cell->file_buffer;
        cell->file_buffer->position_in_map = &grid->buffer_manager->file_map[idx];

        grid->buffer_manager->file_map[idx].prev = NULL;
        grid->buffer_manager->file_map[idx].next = NULL;

        grid->buffer_manager->file_map_tail = grid->buffer_manager->file_map;
    }
    else
    {
        struct file_map *currP = grid->buffer_manager->file_map;
        struct file_map *lastP = grid->buffer_manager->file_map_tail;

        lastP->next = malloc(sizeof(struct file_map));
        if (lastP->next == NULL)
            exit_with_failure("Error in file_buffer_manager: Could not allocate memory for new entry in the file map.\n");

        lastP->next->file_buffer = cell->file_buffer;
        cell->file_buffer->position_in_map = lastP->next;

        lastP->next->prev = lastP;
        lastP->next->next = NULL;
        grid->buffer_manager->file_map_tail = lastP->next;
    }

    grid->buffer_manager->file_map_size++;

    return OK;
}

/* only flushes vectors in metric space */
enum response flush_buffer_to_disk(struct grid *grid, struct cell *cell)
{
    if (cell->file_buffer->buffered_list_size > 0)
    {
        if (cell->filename == NULL)
            exit_with_failure("Error in file_buffer.c: Cannot flush the node to disk. "
                              "It does not have a filename.");

        int full_size = strlen(grid->settings->root_directory) + strlen(cell->filename) + 1;

        char *full_filename = malloc(sizeof(char) * full_size);
        full_filename = strcpy(full_filename, grid->settings->root_directory);
        full_filename = strcat(full_filename, cell->filename);
        full_filename = strcat(full_filename, "\0");

        // if(grid->is_query_grid)
        //     printf("query grid storing in file %s", full_filename);

        // if(grid->is_query_grid == false)
        //     printf("data grid storing in file %s", full_filename);

        COUNT_PARTIAL_OUTPUT_TIME_START
        FILE *vector_file = fopen(full_filename, "a");
        COUNT_PARTIAL_OUTPUT_TIME_END

        if (vector_file == NULL)
        {
            printf("\ncell file name : %s", full_filename);
            exit_with_failure("Error in file_buffer.c: Flushing cell to disk, Could not open the filename.");
        }
        int num_vectors = cell->file_buffer->buffered_list_size;
        int disk_count = cell->file_buffer->disk_count;

        COUNT_PARTIAL_OUTPUT_TIME_START
        for (int i = 0; i < num_vectors; ++i)
        {
            // flush metric vector
            if (!fwrite(cell->file_buffer->mtr_buffered_list[i], sizeof(v_type), grid->settings->mtr_vector_length, vector_file))
                exit_with_failure("Error in file_buffer.c: Could not "
                                  "write metric space vectors to file.\n");
            // flush vector mapping (in pivot space)
            if (!fwrite(cell->file_buffer->ps_buffered_list[i], sizeof(v_type), grid->settings->num_pivots, vector_file))
                exit_with_failure("Error in file_buffer.c: Could not "
                                  "write pivot space vectors to file.\n");
        }
        if (fclose(vector_file))
            exit_with_failure("Error in file_buffer.c: Flushing cell to disk, Could not close file.");
        COUNT_PARTIAL_OUTPUT_TIME_END

        cell->file_buffer->disk_count += num_vectors;

        if (!clear_file_buffer(grid, cell))
            exit_with_failure("Error in file_buffer.c: Flushing cell to disk, Could not clear the buffer.");

        cell->file_buffer->in_disk = true;

        free(full_filename);
    }

    return OK;
}

enum response clear_file_buffer(struct grid *grid, struct cell *cell)
{

    if ((cell->file_buffer) == NULL)
        exit_with_failure("Error in file_buffer.c: Cannot clear a NULL buffer.\n");
    else
    {
        if (cell->file_buffer->mtr_buffered_list != NULL)
            free(cell->file_buffer->mtr_buffered_list);

        if (cell->file_buffer->ps_buffered_list != NULL)
            free(cell->file_buffer->ps_buffered_list);

        cell->file_buffer->mtr_buffered_list = NULL;
        cell->file_buffer->ps_buffered_list = NULL;
        cell->file_buffer->buffered_list_size = 0;
    }

    return OK;
}