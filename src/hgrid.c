#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "../include/hgrid.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_buffer.h"
#include "../include/file_buffer_manager.h"
#include "../include/gsl_matrix.h"
#include "../include/select_pivots.h"

/* initialize grid */
enum response init_grid(const char *root_directory,
                    unsigned int num_pivots,
                    vector *pivots_mtr,
                    vector *pivots_ps,
                    vector *extremity,
                    unsigned int num_levels,
                    unsigned long long num_vectors,
                    unsigned int base,
                    unsigned int mtr_vector_length,
                    double buffered_memory_size,
                    unsigned int max_leaf_size,
                    unsigned int track_vector,
                    struct query_settings * query_settings,
                    struct grid *grid)
{
    // make grid directory
    if (!create_grid_dir(root_directory))
        exit_with_failure("Error in hgrid.c: Couldn't create grid directory!");

    // initialize grid settings
    grid->settings = (struct grid_settings *)malloc(sizeof(struct grid_settings));
    if (grid->settings == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't allocate memory for grid settings!");

    grid->root = NULL;
    grid->first_level = NULL;
    grid->total_records = 0;

    grid->settings->root_directory = root_directory;
    grid->settings->num_pivots = num_pivots; // number of dimensions is equal to the number of pivots
    grid->settings->pivots_mtr = pivots_mtr; // pivot vectors in metric space
    grid->settings->pivots_ps = pivots_ps;   // pivot vectors in pivot space
    grid->settings->pivot_space_extremity = extremity;

    grid->settings->num_levels = num_levels;
    grid->settings->num_leaf_cells = pow(2, num_pivots * num_levels); // 2^(|P| * m) number of cells depends on num_pivots to ensure same length in all edges.
    grid->settings->mtr_vector_length = mtr_vector_length;
    grid->total_records = num_vectors;
    grid->settings->base = base;                              // ex: base = 32 -> one metric vector value is read in 32 bits in a binary file
    grid->settings->max_num_child_cells = pow(2, num_pivots); // equals to number of cells in first level
    grid->settings->vector_size = (base / 8) * mtr_vector_length;

    grid->settings->buffered_memory_size = buffered_memory_size; // amount of memory for file buffers (in MB)
    grid->settings->max_leaf_size = max_leaf_size;               // max number of vectors in one leaf cell
    grid->settings->track_vector = track_vector;
    grid->settings->query_settings = query_settings;

    // track vid
    grid->vid_file = NULL;
    grid->vid_filename = NULL;
    grid->vid_pos_ctr = 0;

    // volume of the pivot space = multiplication of all extremity coordinates.
    grid->settings->pivot_space_volume = 1.0;
    for (int i = 0; i < num_pivots; i++)
    {
        grid->settings->pivot_space_volume *= fabs(extremity->values[i]);
    }

    // leaf cell edge length = (V / num_leaf_cells) ^ 1/|P|
    grid->settings->leaf_cell_edge_length = pow((grid->settings->pivot_space_volume / grid->settings->num_leaf_cells), (1.0 / num_pivots));

    /* 
        Each leaf cell has a file called: 
        level_edge_center(mag,mean)
        02_0.334_(0.33,7.008)
        level: leaf level id (m)
        edge length: edge length of the leaf cell
        mag: magnitude of the center vector of the leaf cell
        number of punctuation marks (underscores:2, parentheses:2, commas:1): total 5
    */

    float edge_length_size = ceil(log10(INT_MAX) + 1);
    float center_vector_size = ceil(log10(INT_MAX)) * 2;
    // float num_vectors_size = ceil(log10(SHRT_MAX) + 1);

    grid->settings->max_filename_size = 2 + edge_length_size + center_vector_size + 5;

    // initialize file buffer manager
    if (!init_file_buffer_manager(grid))
        exit_with_failure("Error in hgrid.c: Could not initialize the file buffer manager for this grid.");

    return OK;
}
/* get extremity vector of the pivot space, pivots must be outliers for extremity to be accurate */
vector *get_extremity(vector *pivot_vectors, unsigned int num_pivots)
{
    // extremity =  farthest vector in the pivot space (holds maximum distance to a pivot)
    vector *pivot_space_extremity = malloc(sizeof(struct vector));
    pivot_space_extremity->values = (v_type *)malloc(sizeof(v_type) * num_pivots);

    v_type max = FLT_MIN;
    for (int i = 0; i < num_pivots; i++)
    {
        for (int j = 0; j < num_pivots; j++)
            if (pivot_vectors[i].values[j] > max)
                max = pivot_vectors[i].values[j];
    }
    for (int i = 0; i < num_pivots; i++)
        pivot_space_extremity->values[i] = max;

    return pivot_space_extremity;
}

/* insert vector in grid */
enum response grid_insert(struct grid *grid, struct inv_index * index, vector *vector)
{
    // get vector mapping in pivot space v -> v'
    struct vector *v_mapping = malloc(sizeof(struct vector));
    if (v_mapping == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't allocate memory for vector mapping.");
    v_mapping->values = malloc(sizeof(v_type) * grid->settings->num_pivots);
    if (v_mapping->values == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't allocate memory for values of vector mapping.");

    map_vector(vector, grid->settings->mtr_vector_length, v_mapping, grid->settings->pivots_mtr, grid->settings->num_pivots);

    // find the closest cell in first level.
    float bsf = FLT_MAX;
    cell *cell = NULL;
    for (int c = 0; c < grid->first_level->num_cells; c++)
    {
        float d = euclidean_distance(v_mapping, grid->first_level->cells[c].center, grid->settings->num_pivots);
        if (d <= bsf)
        {
            bsf = d;
            cell = &grid->first_level->cells[c];
        }
    }

    // loop children of cell untill you find closest leaf cell.
    while (!cell->is_leaf)
    {
        cell = cell_route_to_closest_child(cell, v_mapping, grid->settings->num_pivots);

        if (cell == NULL)
            exit_with_failure("Error in hgrid.c: Could not route to closest child cell.\n");
    }

    // printf("append vector to leaf cell with center :\n");
    // print_vector(cell->center, grid->settings->num_pivots);

    // append vector in metric format
    append_vector_to_cell(grid, index, cell, vector);

    // free memory
    free(v_mapping->values);
    free(v_mapping);

    return OK;
}

/* print grid in console */
void dump_grid_to_console(struct grid *grid)
{
    printf("\n\n\n\t.................................\n");
    printf("\n\n\n\t::          HGRID              ::\n");
    printf("\t.................................\n\n\n");
    printf("|\n|\n|\nv\n");
    level *level = grid->root;
    for (int i = 0; i <= grid->settings->num_levels; i++)
    {
        printf("Level %u:\n", level->id);
        printf("\tNumber of cells = %u\n", level->num_cells);
        printf("\tCell edge length = %f\n\n", level->cell_edge_length);
        for (int j = 0; j < level->num_cells; j++)
        {
            printf("*****************************(Cell %d)********************************\n", j + 1);
            printf("~~Center vector:");
            print_vector(level->cells[j].center, grid->settings->num_pivots);

            if (level->cells[j].children == NULL)
            {
                if (level->is_leaf == false || level->cells[j].is_leaf == false)
                    exit_with_failure("Not leaf level why cell has no children?!");

                printf("(Leaf cell) %s ", level->cells[j].is_leaf ? "true": "false");
                printf("vectors list (total vectors = %u):\n", level->cells[j].cell_size);
                for (int i = 0; i < level->cells[j].cell_size; i++)
                {
                    printf("(%u, %u) \t", level->cells[j].vid[i].table_id, level->cells[j].vid[i].set_pos);
                }
            }
            else
            {
                if (level->is_leaf == true || level->cells[j].is_leaf == true)
                    exit_with_failure("Leaf cell with children?!");
                printf("~~Children:\n");
                for (int k = 0; k < level->cells[j].num_child_cells; k++)
                    print_vector(level->cells[j].children[k].center, grid->settings->num_pivots);
            }
            printf("---\n\n");
        }
        printf("|\n|\n|\nv\n");
        level = level->next;
    }
    printf("\n\t>>>  END OF HGRID  <<<\n\n\n");
}

/* write grid to disk */
enum response grid_write(struct grid *grid)
{
    // (todo) update fct to write root level
    printf(">>> Storing grid : %s\n", grid->settings->root_directory);
    // make root.idx file
    char *root_filename = malloc(sizeof(char) * (strlen(grid->settings->root_directory) + 9));
    if (root_filename == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't allocate memory for root file name.");
    root_filename = strcpy(root_filename, grid->settings->root_directory);
    root_filename = strcat(root_filename, "root.idx\0");

    FILE *root_file = fopen(root_filename, "wb");
    if (root_file == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't open grid file 'root.idx'.");

    free(root_filename);

    // make vid.idx file
    if (grid->settings->track_vector)
    {
        //open vid.edx
        char *vid_filename = malloc(sizeof(char) * (strlen(grid->settings->root_directory) + 8));
        if (vid_filename == NULL)
            exit_with_failure("Error in hgrid.c: Couldn't allocate memory for root.idx");
        strcpy(vid_filename, grid->settings->root_directory);
        strcat(vid_filename, "vid.idx\0");
        grid->vid_file = fopen(vid_filename, "wb");
        if (grid->vid_file == NULL)
            exit_with_failure("Error in hgrid.c: Couldn't open vid.idx");
        grid->vid_pos_ctr = 0;

        free(vid_filename);
    }

    unsigned int num_leaf_cells = grid->settings->num_leaf_cells;
    unsigned int mtr_vector_length = grid->settings->mtr_vector_length;
    unsigned int max_leaf_size = grid->settings->max_leaf_size;
    double buffered_memory_size = grid->settings->buffered_memory_size;
    unsigned long long total_records = grid->total_records;

    // write settings
    fwrite(&num_leaf_cells, sizeof(unsigned long long), 1, root_file);
    fwrite(&buffered_memory_size, sizeof(double), 1, root_file);
    fwrite(&mtr_vector_length, sizeof(unsigned int), 1, root_file);
    fwrite(&max_leaf_size, sizeof(unsigned int), 1, root_file);
    fwrite(&total_records, sizeof(unsigned int), 1, root_file);

    // write levels
    level_write(grid, grid->first_level, root_file);
    fseek(root_file, 0L, SEEK_SET);
    fwrite(&grid->settings->num_levels, sizeof(unsigned int), 1, root_file);

    fclose(root_file);
    fclose(grid->vid_file);
    return OK;
}

/* write level cells to disk */
enum response level_write(struct grid *grid, struct level *level, FILE *file)
{
    fwrite(&(level->is_first), sizeof(unsigned char), 1, file);
    fwrite(&(level->is_leaf), sizeof(unsigned char), 1, file);
    fwrite(&(level->id), sizeof(unsigned int), 1, file);
    fwrite(&(level->num_cells), sizeof(unsigned int), 1, file);
    fwrite(&(level->cell_edge_length), sizeof(unsigned int), 1, file);
    
    if (level->is_leaf)
    {
        for (int c = 0; c < level->num_cells; c++)
        {
            struct cell curr_cell = level->cells[c];
            if (curr_cell.filename != NULL && curr_cell.is_leaf)
            {
                int filename_size = strlen(curr_cell.filename);
                fwrite(&filename_size, sizeof(int), 1, file);
                fwrite(curr_cell.filename, sizeof(char), filename_size, file);
                fwrite(&(curr_cell.cell_size), sizeof(short), 1, file);

                
                if (grid->settings->track_vector)
                {
                    curr_cell.vid_pos = grid->vid_pos_ctr;
                    // save vids to file
                    fwrite(&(curr_cell.vid_pos), sizeof(unsigned int), 1, file);
                    fwrite(curr_cell.vid, sizeof(struct vid), curr_cell.cell_size, grid->vid_file);

                    grid->vid_pos_ctr += curr_cell.cell_size;
                }

                if(curr_cell.file_buffer != NULL)
                    flush_buffer_to_disk(grid, &curr_cell);

                // (todo) update stats
            }
            else
                exit_with_failure("Error in hgrid.c: Cannot write leaf cell to disk without filename!");
        }
    }

    // write next level
    if (!level->is_leaf && level->next != NULL)
        level_write(grid, level->next, file);

    return OK;
}

/* destroy grid */
enum response grid_destroy(struct grid *grid, struct level *level)
{
    //  leaf level
    if (level->is_leaf)
    {
        if (grid->buffer_manager != NULL)
            destroy_buffer_manager(grid);
    }

    // non leaf level
    if (!level->is_leaf)
    {
        grid_destroy(grid, level->next);
    }

    for (int c = level->num_cells - 1; c >= 0; c--)
    {
        // free center
        free(level->cells[c].center->values);

        if (!level->is_first || level->is_root) // center vectors of cells in levels are malloc'd one at a time (check init_levels().)
            free(level->cells[c].center);

        // free filename
        if (level->cells[c].filename != NULL)
            free(level->cells[c].filename);

        // free file buffer
        if (level->cells[c].file_buffer != NULL)
        {
            free(level->cells[c].file_buffer->buffered_list);
            level->cells[c].file_buffer->buffered_list = NULL;
            level->cells[c].file_buffer->buffered_list_size = 0;
            free(level->cells[c].file_buffer);
        }
        // free vid
        if (level->cells[c].vid != NULL)
            free(level->cells[c].vid);
    }

    if (level->is_first) // center vectors of cells in first level are malloc'd all at once (check init_first_level().)
        free(level->cells->center);
    free(level->cells);
    free(level);

    return OK;
}

/* destroy buffer manager */
enum response destroy_buffer_manager(struct grid *grid)
{

    if (grid->buffer_manager != NULL)
    {
        struct file_map *currP;
        struct file_map *temp;

        temp = NULL;
        currP = grid->buffer_manager->file_map;

        while (currP != NULL)
        {
            temp = currP;
            currP = currP->next;
            free(temp);
        }

        free(grid->buffer_manager->memory_array);
        free(grid->buffer_manager);
    }
    return OK;
}
