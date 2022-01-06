#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_buffer.h"
#include "../include/file_buffer_manager.h"
#include "../include/gsl_matrix.h"
#include "../include/select_pivots.h"

/* initialize index */
response init_index(const char *root_directory,
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
                    pexeso_index *index)
{
    // make index directory
    if (!create_index_dir(root_directory))
        exit_with_failure("Error in index.c: Couldn't create index directory!");

    // initialize index settings
    index->settings = (index_settings *)malloc(sizeof(index_settings));
    if (index->settings == NULL)
        exit_with_failure("Error in index.c: Couldn't allocate memory for index settings!");

    index->first_level = NULL;
    index->total_records = 0;

    index->settings->root_directory = root_directory;
    index->settings->num_pivots = num_pivots; // number of dimensions is equal to the number of pivots
    index->settings->pivots_mtr = pivots_mtr; // pivot vectors in metric space
    index->settings->pivots_ps = pivots_ps;   // pivot vectors in pivot space
    index->settings->pivot_space_extremity = extremity;

    index->settings->num_levels = num_levels;
    index->settings->num_leaf_cells = pow(2, num_pivots * num_levels); // 2^(|P| * m) number of cells depends on num_pivots to ensure same length in all edges.
    index->settings->mtr_vector_length = mtr_vector_length;
    index->total_records = num_vectors;
    index->settings->base = base;                              // ex: base = 32 -> one metric vector value is read in 32 bits in a binary file
    index->settings->max_num_child_cells = pow(2, num_pivots); // equals to number of cells in first level
    index->settings->vector_size = (base / 8) * mtr_vector_length;

    index->settings->buffered_memory_size = buffered_memory_size; // amount of memory for file buffers (in MB)
    index->settings->max_leaf_size = max_leaf_size;               // max number of vectors in one leaf cell
    index->settings->track_vector = track_vector;

    // volume of the pivot space = multiplication of all extremity coordinates.
    index->settings->pivot_space_volume = 1.0;
    for (int i = 0; i < num_pivots; i++)
    {
        index->settings->pivot_space_volume *= fabs(extremity->values[i]);
    }

    // leaf cell edge length = (V / num_leaf_cells) ^ 1/|P|
    index->settings->leaf_cell_edge_length = pow((index->settings->pivot_space_volume / index->settings->num_leaf_cells), (1.0 / num_pivots));

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

    index->settings->max_filename_size = 2 + edge_length_size + center_vector_size + 5;

    // initialize file buffer manager
    if (!init_file_buffer_manager(index))
        exit_with_failure("Error in index.c: Could not initialize the file buffer manager for this index.");

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

/* insert vector in index */
response index_insert(pexeso_index *index, vector *vector)
{
    // get vector mapping in pivot space v -> v'
    struct vector *v_mapping = malloc(sizeof(struct vector));
    if (v_mapping == NULL)
        exit_with_failure("Error in index.c: Couldn't allocate memory for vector mapping.");
    v_mapping->values = malloc(sizeof(v_type) * index->settings->num_pivots);
    if (v_mapping->values == NULL)
        exit_with_failure("Error in index.c: Couldn't allocate memory for values of vector mapping.");

    map_vector(vector, index->settings->mtr_vector_length, v_mapping, index->settings->pivots_mtr, index->settings->num_pivots);

    // find the closest cell in first level.
    float bsf = FLT_MAX;
    cell *cell = NULL;
    for (int c = 0; c < index->first_level->num_cells; c++)
    {
        float d = euclidean_distance(v_mapping, index->first_level->cells[c].center, index->settings->num_pivots);
        if (d <= bsf)
        {
            bsf = d;
            cell = &index->first_level->cells[c];
        }
    }

    // loop children of cell untill you find closest leaf cell.
    while (!cell->is_leaf)
    {
        cell = cell_route_to_closest_child(cell, v_mapping, index->settings->num_pivots);

        if (cell == NULL)
            exit_with_failure("Error in index.c: Could not route to closest child cell.\n");
    }

    // printf("append vector to leaf cell with center :\n");
    // print_vector(cell->center, index->settings->num_pivots);

    // append vector in metric format
    append_vector_to_cell(index, cell, vector);

    // free memory
    free(v_mapping->values);
    free(v_mapping);

    return OK;
}

/* print index in console */
void print_index(pexeso_index *index)
{
    printf("DISPLAY PEXESO INDEX...\n");
    printf("|       |       |       \n");
    printf("V       V       V       \n\n\n");
    level *level = index->first_level;
    for (int i = 0; i < index->settings->num_levels; i++)
    {
        printf("Level %u:\n", level->id);
        printf("\tNumber of cells = %u\n", level->num_cells);
        printf("\tCell edge length = %f\n", level->cell_edge_length);
        printf("#####################################################################\n");
        for (int j = 0; j < level->num_cells; j++)
        {
            printf("*****************************(Cell %d)********************************\n", j + 1);
            printf("~~Center vector:");
            print_vector(level->cells[j].center, index->settings->num_pivots);

            if (level->cells[j].children == NULL)
            {
                if (level->is_leaf == false || level->cells[j].is_leaf == false)
                    exit_with_failure("Not leaf level why cell has no children?!");
                printf("(Leaf cell)!!!! ");
                // printf("vectors list (total vectors = %u):\n", level->cells[j].cell_size);
                // for(int i = 0; i < level->cells[j].cell_size; i++)
                // {
                //     printf("(%u, %u) \t", level->cells[j].vid[i].table_id, level->cells[j].vid[i].set_id);
                // }
            }
            else
            {
                if (level->is_leaf == true || level->cells[j].is_leaf == true)
                    exit_with_failure("Leaf cell with children?!");
                printf("~~Children:\n");
                for (int k = 0; k < level->cells[j].num_child_cells; k++)
                    print_vector(level->cells[j].children[k].center, index->settings->num_pivots);
            }
            printf("\nend of cell. ---------------------------------------------------------\n");
        }
        printf("end of level. ########################################################\n");
        level = level->next;
    }
    printf("end of index.\n");
}

/* write index to disk */
enum response index_write(pexeso_index *index)
{
    printf(">>> Storing index : %s\n", index->settings->root_directory);
    // make root.idx file
    char *root_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 9));
    if (root_filename == NULL)
        exit_with_failure("Error in index.c: Couldn't allocate memory for root file name.");
    root_filename = strcpy(root_filename, index->settings->root_directory);
    root_filename = strcat(root_filename, "root.idx\0");

    FILE *root_file = fopen(root_filename, "wb");
    if (root_file == NULL)
        exit_with_failure("Error in index.c: Couldn't open index file 'root.idx'.");

    free(root_filename);

    // make vid.idx file
    if (index->settings->track_vector)
    {
        //open vid.edx
        char *vid_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 8));
        if (vid_filename == NULL)
            exit_with_failure("Error in index.c: Couldn't allocate memory for root.idx");
        strcpy(vid_filename, index->settings->root_directory);
        strcat(vid_filename, "vid.idx\0");
        index->vid_file = fopen(vid_filename, "wb");
        if (index->vid_file == NULL)
            exit_with_failure("Error in index.c: Couldn't open vid.idx");
        index->vid_pos_ctr = 0;

        free(vid_filename);
    }

    unsigned int num_leaf_cells = index->settings->num_leaf_cells;
    unsigned int mtr_vector_length = index->settings->mtr_vector_length;
    unsigned int max_leaf_size = index->settings->max_leaf_size;
    double buffered_memory_size = index->settings->buffered_memory_size;
    unsigned long long total_records = index->total_records;

    // write settings
    fwrite(&num_leaf_cells, sizeof(unsigned long long), 1, root_file);
    fwrite(&buffered_memory_size, sizeof(double), 1, root_file);
    fwrite(&mtr_vector_length, sizeof(unsigned int), 1, root_file);
    fwrite(&max_leaf_size, sizeof(unsigned int), 1, root_file);
    fwrite(&total_records, sizeof(unsigned int), 1, root_file);

    // (todo) write cells and buffers
    // level_write(index, index->first_level, root_file);
    // fseek(root_file, 0L, SEEK_SET);
    // fwrite(&num_leaf_cells, sizeof(unsigned long long), 1, root_file);

    fclose(root_file);

    return OK;
}

/* destroy index */
void index_destroy(struct pexeso_index *index, struct level *level)
{
    //  leaf level
    if (level->is_leaf) 
    {
        if (index->buffer_manager != NULL)
            destroy_buffer_manager(index);
    }

    // non leaf level
    if (!level->is_leaf)
    {
        index_destroy(index, level->next);
    }

    for(int c = level->num_cells - 1; c >= 0; c--)
    {
        // free center vector values
        free(level->cells[c].center->values);
        
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
    // free centers
    free(level->cells->center);
    free(level->cells);
    free(level);
}

/* destroy buffer manager */
enum response destroy_buffer_manager(struct pexeso_index *index)
{

    if (index->buffer_manager != NULL)
    {
        struct file_map *currP;
        struct file_map *temp;

        temp = NULL;
        currP = index->buffer_manager->file_map;

        while (currP != NULL)
        {
            temp = currP;
            currP = currP->next;
            free(temp);
        }

        free(index->buffer_manager->memory_array);
        free(index->buffer_manager);
    }
    return OK;
}