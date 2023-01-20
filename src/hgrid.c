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
#include "../include/stats.h"
#include "../include/inv_index.h"

/* initialize grid */
enum response init_grid(const char *work_dir,
                    unsigned int fft_scale,
                    unsigned int num_pivots,
                    vector *pivots_mtr,
                    vector *pivots_ps,
                    vector *extremity,
                    unsigned int num_levels,
                    unsigned long long num_vectors,
                    unsigned int base,
                    unsigned int mtr_vector_length,
                    double mtr_buffered_memory_size,
                    unsigned int max_leaf_size,
                    unsigned int track_vector,
                    bool is_query_grid,
                    struct query_settings * query_settings,
                    struct grid *grid,
                    unsigned int read_from_disk)
{
    
    // initialize grid settings
    grid->settings = (struct grid_settings *)malloc(sizeof(struct grid_settings));
    if (grid->settings == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't allocate memory for grid settings!");

    grid->root = NULL;
    grid->first_level = NULL;
    grid->total_records = 0;
    grid->is_query_grid  = is_query_grid; 
    
    grid->settings->work_directory = malloc(sizeof(char) * (strlen(work_dir) + 1));
    strcpy(grid->settings->work_directory, work_dir);

    grid->settings->root_directory = malloc(sizeof(char) *( strlen(work_dir) + 8));
    strcpy(grid->settings->root_directory, work_dir);
    if(is_query_grid == true)
        strcat(grid->settings->root_directory, "/Qgrid/\0");
    else
        strcat(grid->settings->root_directory, "/Hgrid/\0");

    grid->settings->fft_scale = fft_scale;
    grid->settings->num_pivots = num_pivots; // number of dimensions is equal to the number of pivots
    grid->settings->pivots_mtr = pivots_mtr; // pivot vectors in metric space
    grid->settings->pivots_ps = pivots_ps;   // pivot vectors in pivot space
    grid->settings->pivot_space_extremity = extremity;

    grid->settings->num_levels = num_levels;
    grid->settings->num_leaf_cells = pow(2, num_pivots * num_levels); // 2^(|P| * m) number of cells depends on num_pivots to ensure same length in all edges.
    grid->settings->mtr_vector_length = mtr_vector_length;
    grid->settings->base = base;                              // ex: base = 32 -> one metric vector value is read in 32 bits in a binary file
    grid->settings->max_num_child_cells = pow(2, num_pivots); // equals to number of cells in first level
    grid->settings->vector_size = (base / 8) * mtr_vector_length;

    grid->settings->mtr_buffered_memory_size = mtr_buffered_memory_size; // amount of memory for file buffers (metric vectors) (in MB)

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
        number of punctuation marks (underscores:2, parentheses:2): total 4
        commas: total (num_pivots - 1)
    */

    float edge_length_size = ceil(log10(INT_MAX) + 1);
    float center_vector_size = ceil(log10(INT_MAX)) * grid->settings->num_pivots;
    // float num_vectors_size = ceil(log10(SHRT_MAX) + 1);

    grid->settings->max_filename_size = 2 + edge_length_size + center_vector_size + 4 + grid->settings->num_pivots - 1;

    // make grid directory
    if(read_from_disk == 0 && !is_query_grid)
        if (!create_grid_dir(grid->settings->root_directory))
            exit_with_failure("Error in hgrid.c: Couldn't create grid directory!");

    // initialize file buffer manager
    if (!init_file_buffer_manager(grid))
        exit_with_failure("Error in hgrid.c: Could not initialize the file buffer manager for this grid.");

    return OK;
}


/* init statistics */
enum response init_grid_stats(struct grid * grid)
{
    grid->stats = malloc(sizeof(struct stats_info));
    if (grid->stats == NULL)
        exit_with_failure("Error in hgrid.c: Could not allocate memory for stats structure.");
    
    grid->stats->total_cells_count = 0;
    grid->stats->leaf_cells_count = 0;
    grid->stats->empty_leaf_cells_count = 0;

    grid->stats->loaded_files_count = 0;
    grid->stats->loaded_query_files_count = 0;

    grid->stats->loaded_sets_count = 0;
    grid->stats->loaded_query_sets_count = 0; 

    grid->stats->loaded_files_size = 0;
    grid->stats->loaded_query_files_size = 0;
    grid->stats->loaded_vec_count = 0;
    grid->stats->loaded_qvec_count = 0;

    grid->stats->out_of_ps_space_vec_count = 0;
    grid->stats->out_of_ps_space_qvec_count = 0;

    grid->stats->total_queries_count = 0;

    for(int i = 0; i < 7; i++)
        grid->stats->used_lemmas_count[i] = 0;

    grid->stats->filtered_cells_count = 0;
    grid->stats->visited_cells_count = 0;
    grid->stats->visited_matching_cells_count = 0;
    grid->stats->visited_candidate_cells_count = 0;

    grid->stats->filtered_vectors_count = 0;
    grid->stats->matching_vectors_count = 0;
    grid->stats->checked_vectors_in_ps_count = 0;
    grid->stats->checked_vectors_in_mtr_count = 0;

    // timers 
    grid->stats->idx_append_vec_to_leaf_total_time = 0;	
	grid->stats->idx_append_vec_to_leaf_input_time = 0;
	grid->stats->idx_append_vec_to_leaf_output_time = 0;
	grid->stats->idx_append_vec_to_leaf_cpu_time = 0;

    grid->stats->total_input_time = 0;
    grid->stats->total_output_time = 0;
    grid->stats->total_parse_time = 0;
    grid->stats->total_query_time = 0;
    grid->stats->total_pivot_selection_time = 0;
    grid->stats->total_time = 0;

    grid->stats->grid_building_total_time = 0;
    grid->stats->grid_writing_total_time = 0;

    return OK;  
}


/* init statistics */
enum response reset_grid_stats(struct grid * grid)
{
    grid->stats->total_cells_count = 0;
    grid->stats->leaf_cells_count = 0;
    grid->stats->empty_leaf_cells_count = 0;

    grid->stats->loaded_files_count = 0;
    grid->stats->loaded_query_files_count = 0;

    grid->stats->loaded_sets_count = 0;
    grid->stats->loaded_query_sets_count = 0; 

    grid->stats->loaded_files_size = 0;
    grid->stats->loaded_query_files_size = 0;
    grid->stats->loaded_vec_count = 0;
    grid->stats->loaded_qvec_count = 0;

    grid->stats->out_of_ps_space_qvec_count = 0;

    for(int i = 0; i < 7; i++)
        grid->stats->used_lemmas_count[i] = 0;

    grid->stats->filtered_cells_count = 0;
    grid->stats->visited_cells_count = 0;
    grid->stats->visited_matching_cells_count = 0;
    grid->stats->visited_candidate_cells_count = 0;

    grid->stats->filtered_vectors_count = 0;
    grid->stats->matching_vectors_count = 0;
    grid->stats->checked_vectors_in_ps_count = 0;
    grid->stats->checked_vectors_in_mtr_count = 0;

    // timers 
    grid->stats->idx_append_vec_to_leaf_total_time = 0;	
	grid->stats->idx_append_vec_to_leaf_input_time = 0;
	grid->stats->idx_append_vec_to_leaf_output_time = 0;
	grid->stats->idx_append_vec_to_leaf_cpu_time = 0;

    grid->stats->total_input_time = 0;
    grid->stats->total_output_time = 0;
    grid->stats->total_parse_time = 0;
    grid->stats->total_query_time = 0;

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
    
    v_mapping->set_id = vector->set_id;
    v_mapping->table_id = vector->table_id;
    v_mapping->pos = vector->pos;

    v_mapping->values = malloc(sizeof(v_type) * grid->settings->num_pivots);
    if (v_mapping->values == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't allocate memory for values of vector mapping.");

    map_vector(vector, grid->settings->mtr_vector_length, v_mapping, grid->settings, grid->is_query_grid);
    
    
    // find the closest cell in first level.
    float bsf = FLT_MAX;
    cell *cell = NULL;
    for (int c = 0; c < grid->first_level->num_cells; c++)
    {
        float d = euclidean_distance_cmp(v_mapping, grid->first_level->cells[c].center, grid->settings->num_pivots);
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

    // printf("append vector to leaf cell with mapping :\n");
    // print_vector(v_mapping, grid->settings->num_pivots);

    // append vector in metric format
    if(!append_vector_to_cell(grid, index, cell, vector, v_mapping))
        exit_with_failure("Error in hgrid.c: couldn't append vector to cell.");

    // if(grid->is_query_grid)
    //     print_vector(v_mapping, grid->settings->num_pivots);
    // free memory
    free(v_mapping->values);
    free(v_mapping);

    return OK;
}

/* print grid in console */
void dump_grid_to_console(struct grid *grid)
{
    printf("\n\n\n\t.................................\n");
    printf("\t::          GRID              ::\n");
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
                    printf("(t %u, c %u, pos %u, size %u) \t", level->cells[j].vid[i].table_id, level->cells[j].vid[i].set_id, level->cells[j].vid[i].pos, level->cells[j].vid[i].set_size);
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
    printf("\n\t>>>  END OF GRID  <<<\n\n\n");
}

/* write grid to disk */
enum response grid_write(struct grid *grid)
{
    // (todo) update fct to write root level
    printf("Storing grid : %s\n", grid->settings->root_directory);
    // make root.idx file
    char *root_filename = malloc(sizeof(char) * (strlen(grid->settings->root_directory) + 9));
    if (root_filename == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't allocate memory for root file name.");
    root_filename = strcpy(root_filename, grid->settings->root_directory);
    root_filename = strcat(root_filename, "root.idx\0");

    COUNT_PARTIAL_OUTPUT_TIME_START
    FILE *root_file = fopen(root_filename, "wb");
    COUNT_PARTIAL_OUTPUT_TIME_END
    if (root_file == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't open grid file 'root.idx'.");

    // make vid.idx file
    if (grid->settings->track_vector)
    {
        //open vid.edx
        char *vid_filename = malloc(sizeof(char) * (strlen(grid->settings->root_directory) + 8));
        if (vid_filename == NULL)
            exit_with_failure("Error in hgrid.c: Couldn't allocate memory for root.idx");
        strcpy(vid_filename, grid->settings->root_directory);
        strcat(vid_filename, "vid.idx\0");
        COUNT_PARTIAL_OUTPUT_TIME_START
        grid->vid_file = fopen(vid_filename, "wb");
        COUNT_PARTIAL_OUTPUT_TIME_END
        if (grid->vid_file == NULL)
            exit_with_failure("Error in hgrid.c: Couldn't open vid.idx");
        grid->vid_pos_ctr = 0;

        free(vid_filename);
    }

    unsigned int num_levels = grid->settings->num_levels;
    unsigned int num_pivots = grid->settings->num_pivots;
    unsigned int mtr_vector_length = grid->settings->mtr_vector_length;
    v_type *pivot_space_extremity = grid->settings->pivot_space_extremity->values;
    struct vector *pivots_mtr = grid->settings->pivots_mtr;
    struct vector *pivots_ps = grid->settings->pivots_ps;

    unsigned int num_leaf_cells = grid->settings->num_leaf_cells;
    unsigned int max_leaf_size = grid->settings->max_leaf_size;
    double mtr_buffered_memory_size = grid->settings->mtr_buffered_memory_size;
    double ps_buffered_memory_size = grid->settings->ps_buffered_memory_size;
    unsigned long long total_records = grid->total_records;

    // write settings
    COUNT_PARTIAL_OUTPUT_TIME_START
    fwrite(&num_levels, sizeof(unsigned int), 1, root_file);
    fwrite(&num_pivots, sizeof(unsigned int), 1, root_file);
    fwrite(&mtr_vector_length, sizeof(unsigned int), 1, root_file);
    fwrite(pivot_space_extremity, sizeof(v_type), num_pivots, root_file); // write extrimity
    for(int p = 0; p < num_pivots; p++) // write pivots in metric space
    {
        fwrite(pivots_mtr[p].values, sizeof(v_type), mtr_vector_length, root_file);
    }
    for(int p = 0; p < num_pivots; p++) // write pivots in pivot space
    {
        fwrite(pivots_ps[p].values, sizeof(v_type), num_pivots, root_file);
    }
    fwrite(&num_leaf_cells, sizeof(unsigned int), 1, root_file);
    fwrite(&mtr_buffered_memory_size, sizeof(double), 1, root_file);
    fwrite(&ps_buffered_memory_size, sizeof(double), 1, root_file);
    fwrite(&max_leaf_size, sizeof(unsigned int), 1, root_file);
    fwrite(&total_records, sizeof(unsigned long long), 1, root_file);
    COUNT_PARTIAL_OUTPUT_TIME_END

    // write levels
    struct level * leaf_level, *current_level = grid->first_level;
    while(!current_level->is_leaf)
        current_level = current_level->next;
    leaf_level = current_level;
    write_leaf_level(grid, leaf_level, root_file);

    COUNT_PARTIAL_OUTPUT_TIME_START
    fseek(root_file, 0L, SEEK_SET);
    fwrite(&grid->settings->num_levels, sizeof(unsigned int), 1, root_file);
    COUNT_PARTIAL_OUTPUT_TIME_END


    free(root_filename);
    COUNT_PARTIAL_OUTPUT_TIME_START
    fclose(root_file);
    fclose(grid->vid_file);
    COUNT_PARTIAL_OUTPUT_TIME_END
    return OK;
}

/* read grid from disk */
struct grid * grid_read(const char * grid_directory, struct query_settings * query_settings)
{
    struct grid * grid = (struct grid *) malloc(sizeof(struct grid));
    if (grid == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for grid!");
    
    printf("Reading grid from : %s\n", grid_directory);
    // read root.idx file
    char *root_filename = malloc(sizeof(char) * (strlen(grid_directory) + 9));
    if (root_filename == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't allocate memory for root file name.");
    root_filename = strcpy(root_filename, grid_directory);
    root_filename = strcat(root_filename, "root.idx\0");

    COUNT_PARTIAL_INPUT_TIME_START
    FILE *root_file = fopen(root_filename, "rb");
    COUNT_PARTIAL_INPUT_TIME_END
    if (root_file == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't open grid file 'root.idx'.");

    // read vid.idx file
    char *vid_filename = malloc(sizeof(char) * (strlen(grid_directory) + 8));
        if (vid_filename == NULL)
            exit_with_failure("Error in hgrid.c: Couldn't allocate memory for root.idx");
    strcpy(vid_filename,grid_directory);
    strcat(vid_filename, "vid.idx\0");

    COUNT_PARTIAL_INPUT_TIME_START
    FILE *vid_file = fopen(vid_filename, "rb");
    COUNT_PARTIAL_INPUT_TIME_END
    if (vid_file == NULL)
        exit_with_failure("Error in hgrid.c: Couldn't open vid.idx");
    
    unsigned int num_levels = 0;
    unsigned int num_pivots = 0;
    unsigned int mtr_vector_length = 0;
    struct vector *pivot_space_extremity = malloc(sizeof(struct vector));
    struct vector *pivots_mtr = NULL;
    struct vector *pivots_ps = NULL;

    unsigned int num_leaf_cells = 0;
    unsigned int max_leaf_size = 0;
    double mtr_buffered_memory_size = 0;
    double ps_buffered_memory_size = 0;
    unsigned long long total_records = 0;


    // read settings
    printf("\n\nRead grid settings... ");
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(root_file, 0L, SEEK_SET);    
    fread(&num_levels, sizeof(unsigned int), 1, root_file);
    fread(&num_pivots, sizeof(unsigned int), 1, root_file);
    fread(&mtr_vector_length, sizeof(unsigned int), 1, root_file);

    pivot_space_extremity->values = malloc(sizeof(v_type) * num_pivots); 
    if(pivot_space_extremity->values == NULL) exit_with_failure("Error in hgrid.c:  Could not allocate memory for extrimity.");
    fread(pivot_space_extremity->values, sizeof(v_type), num_pivots, root_file); // read extrimity

    pivots_mtr = malloc(sizeof(struct vector) * num_pivots);
    pivots_ps = malloc(sizeof(struct vector) * num_pivots);
    if(pivots_mtr == NULL || pivots_ps == NULL) exit_with_failure("Error in hgrid.c:  Could not allocate memory for pivot vectors.");

    for(int p = 0; p < num_pivots; p++) // write pivots in metric space
    {
        pivots_mtr[p].values = malloc(sizeof(v_type) * mtr_vector_length);
        if(pivots_mtr[p].values == NULL) exit_with_failure("Error in hgrid.c:  Could not allocate memory for pivot values.");
        fread(pivots_mtr[p].values, sizeof(v_type), mtr_vector_length, root_file);
    }
    for(int p = 0; p < num_pivots; p++) // write pivots in pivot space
    {
        pivots_ps[p].values = malloc(sizeof(v_type) * num_pivots);
        if(pivots_ps[p].values == NULL) exit_with_failure("Error in hgrid.c:  Could not allocate memory for pivot values.");
        fread(pivots_ps[p].values, sizeof(v_type), num_pivots, root_file);
    }
    fread(&num_leaf_cells, sizeof(unsigned int), 1, root_file);
    fread(&mtr_buffered_memory_size, sizeof(double), 1, root_file);
    fread(&ps_buffered_memory_size, sizeof(double), 1, root_file);
    fread(&max_leaf_size, sizeof(unsigned int), 1, root_file);
    fread(&total_records, sizeof(unsigned long long), 1, root_file);
    COUNT_PARTIAL_INPUT_TIME_END

    printf("\n\nInit grid... ");
    if (!init_grid(grid_directory, 0, num_pivots, pivots_mtr, pivots_ps, pivot_space_extremity, 
                    num_levels, total_records, 32, mtr_vector_length, 
                    mtr_buffered_memory_size, max_leaf_size, 1, 
                    false, query_settings, grid, 1))
        exit_with_failure("Error in main.c: Couldn't initialize grid!");
    
    printf("\n\nInit grid stats... ");
    if(!init_grid_stats(grid))
        exit_with_failure("Error in main.c: Couldn't initialize grid stats!");
    printf("(OK)\n");


    // Build levels
    printf("\n\nRebuilding levels... ");
    if(!init_root(grid))
        exit_with_failure("Error in main.c: Couldn't initialize root level!");

    if (!init_first_level(grid))
        exit_with_failure("Error in main.c: Couldn't initialize first level!");

    if (!init_levels(grid))
        exit_with_failure("Error in main.c: Couldn't initialize grid levels!");
    printf("(OK)\n");

    printf("\n\nReading leaf cells... ");
    struct level * leaf_level, *current_level = grid->first_level;
    while(!current_level->is_leaf)
        current_level = current_level->next;
    leaf_level = current_level;
    printf("(OK)\n");

    // read leaf level
    read_leaf_level(grid, root_file, leaf_level);
    
    if(!count_empty_leaf_cells(grid))
            exit_with_failure("Error in main.c: Couldn't count NB of empty leaf cells.");


    // print grid settings
    print_grid_settings(grid);

    COUNT_PARTIAL_INPUT_TIME_START
    fclose(root_file);
    fclose(vid_file);
    free(root_filename);
    free(vid_filename);
    COUNT_PARTIAL_INPUT_TIME_END

    return grid;
}

/* read level cells from disk */
enum response read_leaf_level(struct grid *grid, FILE *root_file, struct level * level)
{
    bool is_root = true;
    bool is_first = true;
    bool is_leaf = true;
    unsigned int id = 0, num_cells = 0;
    float cell_edge_length = 0.0;

    COUNT_PARTIAL_INPUT_TIME_START
    fread(&is_root, sizeof(unsigned char), 1, root_file);
    fread(&is_first, sizeof(unsigned char), 1, root_file);
    fread(&is_leaf, sizeof(unsigned char), 1, root_file);
    fread(&id, sizeof(unsigned int), 1, root_file);
    fread(&num_cells, sizeof(unsigned int), 1, root_file);
    fread(&cell_edge_length, sizeof(unsigned int), 1, root_file);
    COUNT_PARTIAL_INPUT_TIME_END
    

    int filename_size = 0;
    char filename[255] = "";
    unsigned int vid_pos = 0;
    struct vector temp_center;
    struct vid * vid = NULL;
    temp_center.values = malloc(sizeof(v_type) * grid->settings->num_pivots);
    if(temp_center.values == NULL)
        exit_with_failure("Error in hgrid.c: Couldnt allocate memory for temp vector.");
        
    for (int c = 0; c < level->num_cells; c++)
    {
        struct cell * cell = &(level->cells[c]);
        COUNT_PARTIAL_INPUT_TIME_START
        fread(&filename_size, sizeof(int), 1, root_file);
        fread(filename, sizeof(char), filename_size, root_file);
        fread(temp_center.values, sizeof(v_type), grid->settings->num_pivots, root_file);
        fread(&(cell->cell_size), sizeof(unsigned int), 1, root_file);

        if (grid->settings->track_vector)
        {
            if(cell->cell_size != 0)
            {
                cell->vid = realloc(cell->vid, sizeof(struct vid) * cell->cell_size);
                if(cell->vid == NULL)
                    exit_with_failure("Error in hgrid.c: Couldnt allocate memory for cell vid vector.");
            
                fread(&(cell->vid_pos), sizeof(unsigned int), 1, root_file);
                fread(cell->vid, sizeof(struct vid), cell->cell_size, root_file);

                // printf("<> cell %d, size %d:\n", cell->id, cell->cell_size);
                // printf("+ vector ids:\n");
                // for(int v = 0; v < cell->cell_size; v++)
                // {
                //     printf("(%u, %u, %u)\t", cell->vid[v].table_id, cell->vid[v].set_id, cell->vid[v].pos);
                // }
                // printf("\n");
            }
            if (!file_buffer_init(cell))
            exit_with_failure("Error in hgrid.c: Could not initialize the file buffer \
                           for the node.\n");

            if(cell->file_buffer == NULL)
                exit_with_failure("Error in hgrid.c: Could not initialize the file buffer \
                            for the node.\n");

            cell->file_buffer->in_disk = true;
            cell->file_buffer->disk_count = cell->cell_size;
            cell->file_buffer->buffered_list_size = 0; 
        }
        COUNT_PARTIAL_INPUT_TIME_END
    }
    free(temp_center.values);
    return OK;
}

/* write level cells to disk */
enum response write_leaf_level(struct grid *grid, struct level *level, FILE *root_file)
{
    COUNT_PARTIAL_OUTPUT_TIME_START
    fwrite(&(level->is_root), sizeof(unsigned char), 1, root_file);
    fwrite(&(level->is_first), sizeof(unsigned char), 1, root_file);
    fwrite(&(level->is_leaf), sizeof(unsigned char), 1, root_file);
    fwrite(&(level->id), sizeof(unsigned int), 1, root_file);
    fwrite(&(level->num_cells), sizeof(unsigned int), 1, root_file);
    fwrite(&(level->cell_edge_length), sizeof(unsigned int), 1, root_file);
    COUNT_PARTIAL_OUTPUT_TIME_END

    if (!level->is_leaf)
        exit_with_failure("Error in hgrid.c: Cannot write non leaf level!");

    // if (level->is_leaf)   
    // {
        for (int c = 0; c < level->num_cells; c++)
        {
            
            struct cell curr_cell = level->cells[c];
            
            if (curr_cell.filename != NULL && curr_cell.is_leaf)
            {
                int filename_size = strlen(curr_cell.filename);
                COUNT_PARTIAL_OUTPUT_TIME_START
                fwrite(&filename_size, sizeof(int), 1, root_file);
                fwrite(curr_cell.filename, sizeof(char), filename_size, root_file);
                fwrite(curr_cell.center->values, sizeof(v_type), grid->settings->num_pivots, root_file);
                fwrite(&(curr_cell.cell_size), sizeof(unsigned int), 1, root_file);
                COUNT_PARTIAL_OUTPUT_TIME_END

                // (todo) flush_buffer_to_disk
                if (grid->settings->track_vector)
                {
                    curr_cell.vid_pos = grid->vid_pos_ctr;
                    // save vids to file
                    COUNT_PARTIAL_OUTPUT_TIME_START
                    if(curr_cell.cell_size != 0)
                    {
                        fwrite(&(curr_cell.vid_pos), sizeof(unsigned int), 1, root_file);
                        fwrite(curr_cell.vid, sizeof(struct vid), curr_cell.cell_size, root_file);

                        // printf("<> cell %d, size %d:\n", curr_cell.id, curr_cell.cell_size);
                        // printf("+ vector ids:\n");
                        // for(int v = 0; v < curr_cell.cell_size; v++)
                        // {
                        //     printf("(%u, %u, %u)\t", curr_cell.vid[v].table_id, curr_cell.vid[v].set_id, curr_cell.vid[v].pos);
                        // }
                        // printf("\n");
                    }
                    COUNT_PARTIAL_OUTPUT_TIME_END

                    grid->vid_pos_ctr += curr_cell.cell_size;
                }

                if(curr_cell.file_buffer != NULL)
                    flush_buffer_to_disk(grid, &curr_cell);
                // (todo) update stats
            }
            else
                exit_with_failure("Error in hgrid.c: Cannot write leaf cell to disk without filename!");
        }
    // }

    // write next level
    // if (!level->is_leaf && level->next != NULL)
    //     write_leaf_level(grid, level->next, root_file);

    return OK;
}

/* destroy grid levels*/
enum response grid_destroy_level(struct grid *grid, struct level *level)
{
    //  leaf level
    if (level->is_leaf == true)
    {
        if (grid->buffer_manager != NULL)
            destroy_buffer_manager(grid);
    }

    // non leaf level
    if (!level->is_leaf)
    {
        grid_destroy_level(grid, level->next);
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
            free(level->cells[c].file_buffer->mtr_buffered_list);
            free(level->cells[c].file_buffer->ps_buffered_list);
            level->cells[c].file_buffer->mtr_buffered_list = NULL;
            level->cells[c].file_buffer->ps_buffered_list = NULL;
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

/* destroy grid */
enum response grid_destroy(struct grid *grid)
{
    if(grid->is_query_grid)
        exit_with_failure("Error in hgrid.c: Calling grid destroy for query grid, please call query_grid_destroy() instead.");
    
    grid_destroy_level(grid, grid->root);

    for(int p = grid->settings->num_pivots - 1; p >= 0; p--)
        free(grid->settings->pivots_mtr[p].values);
    free(grid->settings->pivots_mtr);

    for(int p = grid->settings->num_pivots - 1; p >= 0; p--)
        free(grid->settings->pivots_ps[p].values);
    free(grid->settings->pivots_ps);
    
    free(grid->settings->pivot_space_extremity->values);
    free(grid->settings->pivot_space_extremity);
    free(grid->settings->query_settings);

    free(grid->settings->root_directory);
    free(grid->settings->work_directory);
    // destroy settings
    free(grid->settings);
    free(grid->stats);
    free(grid);

    return OK;
}

/* destroy query grid */
enum response query_grid_destroy(struct grid *grid)
{
    // query grid shares settings with data grid, only free levels and the grid
    grid_destroy_level(grid, grid->root);

    free(grid->settings->root_directory);
    free(grid->settings->work_directory);

    free(grid->settings);
    free(grid);

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

        free(grid->buffer_manager->mtr_memory_array);
        free(grid->buffer_manager->ps_memory_array);
        free(grid->buffer_manager);
        grid->buffer_manager = NULL;
    }
    return OK;
}

/* print stats */
void print_grid_stats(struct grid * grid)
{
    printf("\n\n\n(d) Datalake:\t-------------------------------------\n\n\n");
    printf("Loaded_files_count\t%ld\n", 
        grid->stats->loaded_files_count);
    
    printf("Loaded_sets_count\t%ld\n", 
        grid->stats->loaded_sets_count);

    printf("Loaded_vec_count\t%ld\n", 
        grid->stats->loaded_vec_count);
    
    printf("Loaded_files_size\t%.3f (in GB)\n", 
        grid->stats->loaded_files_size / (1024 * 1024 * 1024));
    
    printf("\n\n\n(q) Query:\t-------------------------------------\n\n\n");

    printf("Loaded_query_files_count\t%ld\n", 
        grid->stats->loaded_query_files_count);

    printf("Loaded_query_sets_count\t%ld\n", 
        grid->stats->loaded_query_sets_count);
    
    printf("Loaded_qvec_count\t%ld\n", 
        grid->stats->loaded_qvec_count);
    
    printf("Loaded_query_files_size\t%.3f (in GB)\n", 
        grid->stats->loaded_query_files_size / (1024 * 1024 * 1024));

    printf("\n\n\n(s) Settings:\t-------------------------------------\n\n\n");

    printf("Pivots_count\t%d\n", 
        grid->settings->num_pivots);
    
    printf("Levels_count\t%d\n", 
        grid->settings->num_levels);

    printf("Leaf_cells_count\t%d\n", 
        grid->settings->num_leaf_cells);
    
    printf("Empty_leaf_cells_count\t%d\n", 
        grid->stats->empty_leaf_cells_count);
    
    printf("\n\n\n(!) Warnings:\t-------------------------------------\n\n\n");
    printf("Out_of_ps_space_vec_count\t%ld\n", 
        grid->stats->out_of_ps_space_vec_count);

    printf("Out_of_ps_space_qvec_count\t%ld\n", 
        grid->stats->out_of_ps_space_qvec_count);

    printf("\n\n\n(t) Time measures in seconds:\t-------------------------------------\n\n\n");

    printf("Total_time\t%lf\n",
         grid->stats->total_time / 1000000);

    printf("Grid_building_total_time\t%lf\n",
         grid->stats->grid_building_total_time / 1000000);
        
    printf("Grid_building_input_time\t%lf\n",
         grid->stats->grid_building_input_time / 1000000);

    printf("Grid_building_output_time\t%lf\n",
         grid->stats->grid_building_output_time / 1000000);

    printf("Total_pivot_selection_time\t%lf\n", 
        grid->stats->total_pivot_selection_time / 1000000);

    printf("Total_query_time\t%lf\n", 
        grid->stats->total_query_time / 1000000);

        printf("\n\n\n(!) Lemma usage:\t-------------------------------------\n\n\n");
    for(int i = 0; i < 7; i++)
        printf("Lemma %d:\t%u\n", i+1, grid->stats->used_lemmas_count[i]);

    printf("\n\n\n(!) Data access:\t-------------------------------------\n\n\n");
    printf("Filtered_cells_count\t%d\n", grid->stats->filtered_cells_count);
    printf("Visited_cells_count\t%d\n", grid->stats->visited_cells_count);
    printf("Visited_matching_cells_count\t%d\n", grid->stats->visited_matching_cells_count);
    printf("Visited_candidate_cells_count\t%d\n\n", grid->stats->visited_candidate_cells_count);

    printf("filtered_vectors_count\t%ld\n", grid->stats->filtered_vectors_count);
    printf("matching_vectors_count\t%ld\n", grid->stats->matching_vectors_count);
    printf("Checked_vectors_in_ps_count\t%ld\n", grid->stats->checked_vectors_in_ps_count);
    printf("Checked_vectors_in_mtr_count\t%ld\n\n", grid->stats->checked_vectors_in_mtr_count);

    printf("Count_add_mpair\t%d\n", grid->stats->count_add_mpair);
    printf("Count_add_cpair\t%d\n", grid->stats->count_add_cpair);
    printf("Count_dist_calc\t%ld\n", grid->stats->count_dist_calc);


    printf("fft_scale,\ttau,\tT,\t%%filled_cells,\t#vec_outps,\t#filtered_vectors,\t#matching_vectors"
            ",\t#checked_pvectors,\t#checked_mvectors,\t#visited_cells,\t#visited_mcells,\t#visited_ccells,\t#filtered_cells\n");
    printf("%d,\t%f,\t%f,\t%.2f%%,\t%ld,\t%ld,\t%ld,\t%ld,\t%ld"
            ",\t%d,\t%d,\t%d,\t%d\n", grid->settings->fft_scale, 
            grid->settings->query_settings->dist_threshold/2, 
            grid->settings->query_settings->join_threshold,
            100.0 * (((float)grid->settings->num_leaf_cells - (float)grid->stats->empty_leaf_cells_count)/(float)grid->settings->num_leaf_cells),
            grid->stats->out_of_ps_space_vec_count,
            grid->stats->filtered_vectors_count, grid->stats->matching_vectors_count,
            grid->stats->checked_vectors_in_ps_count, grid->stats->checked_vectors_in_mtr_count,
            grid->stats->visited_cells_count, grid->stats->visited_matching_cells_count, grid->stats->visited_candidate_cells_count,
            grid->stats->filtered_cells_count
            );
    printf("\n\n\n"); 
}

/* print settings */
void print_grid_settings(struct grid * grid)
{
    printf("\n\n--- --- HGRID SETTINGS --- ---\n\n");
    printf("m:\t\t%d\n", grid->settings->num_levels);
    printf("P:\t\t%d\n", grid->settings->num_pivots);
    printf("mtr_dim:\t\t%d\n", grid->settings->mtr_vector_length);
    printf("Extrimity:\t\t\n");
    print_vector(grid->settings->pivot_space_extremity, grid->settings->num_pivots);
    printf("Pivots in mtr space:\t\t\n");
    for(int p = 0; p < grid->settings->num_pivots; p++)
    {
        printf("P%d:\n", p+1);
        print_vector(&(grid->settings->pivots_mtr[p]), grid->settings->mtr_vector_length);
    }
    printf("Pivots in pivot space:\t\t\n");
    for(int p = 0; p < grid->settings->num_pivots; p++)
    {
        printf("P%d:\n", p+1);
        print_vector(&(grid->settings->pivots_ps[p]), grid->settings->num_pivots);
    }
    printf("#Leaf_cells:\t\t%d\n", grid->settings->num_leaf_cells);
    printf("#Empty_leaf_cells:\t\t%d\n", grid->stats->empty_leaf_cells_count);
    printf("#Max_leaf_size:\t\t%d\n", grid->settings->max_leaf_size);
    printf("Mtr_buffered_memory_size:\t\t%.3f\n", grid->settings->mtr_buffered_memory_size);
    printf("Ps_buffered_memory_size:\t\t%.3f\n", grid->settings->ps_buffered_memory_size);
    printf("Total_records:\t\t%lld\n", grid->total_records);
    printf("\n\n--- --- --- --- --- --- ---\n\n");
}

/* count empty leaf cells */
enum response count_empty_leaf_cells(struct grid * grid)
{
    int num_empty_leaves = 0;

    struct level * level = grid->root;
    while(!level->is_leaf) // go to leaf level
        level = level->next;

    for(int c = 0; c < level->num_cells; c++)
    {
        if(level->cells[c].cell_size == 0)
            num_empty_leaves += 1;

        else if(level->cells[c].cell_size < 0)
            exit_with_failure("Error in hgrid.c: counting empty leafs, found leaf with negative size!");
    }

    grid->stats->empty_leaf_cells_count = num_empty_leaves;
    return num_empty_leaves;
}