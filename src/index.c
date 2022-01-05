#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_buffer.h"
#include "../include/file_buffer_manager.h"
#include "../include/gsl_matrix.h"
#include "../include/select_pivots.h"

/* initialize index */
response init_index(const char * root_directory,
                unsigned int num_pivots,
                vector * extremity,
                unsigned int num_levels,
                unsigned long long num_vectors,
                unsigned int base,
                unsigned int mtr_vector_length,
                double buffered_memory_size,
                unsigned int max_leaf_size,
                pexeso_index * index)
{
    // initialize index settings
    index->settings = (index_settings *)malloc(sizeof(index_settings));
    if (index->settings == NULL)
        exit_with_failure("Error in index.c: Couldn't allocate memory for index settings!");


    index->first_level = NULL;
    index->total_records = 0;

    index->settings->root_directory = root_directory;
    index->settings->num_pivots = num_pivots; // number of dimensions is equal to the number of pivots
    index->settings->pivot_space_extremity = extremity;
    index->settings->num_levels = num_levels;
    index->settings->num_leaf_cells = pow(2, num_pivots * num_levels); // 2^(|P| * m) number of cells depends on num_pivots to ensure same length in all edges.
    index->settings->mtr_vector_length = mtr_vector_length;
    index->total_records = num_vectors;
    index->settings->base = base; // ex: base = 32 -> one metric vector value is read in 32 bits in a binary file 
    index->settings->max_num_child_cells = pow(2, num_pivots); // equals to number of cells in first level
    index->settings->vector_size = (base/8) * mtr_vector_length;
    
    index->settings->buffered_memory_size = buffered_memory_size; // amount of memory for file buffers (in MB)
    index->settings->max_leaf_size = max_leaf_size; // max number of vectors in one leaf cell

    // volume of the pivot space = multiplication of all extremity coordinates.
    index->settings->pivot_space_volume = 1.0;
    for(int i = 0; i < num_pivots; i++)
    {
        index->settings->pivot_space_volume *= fabs(extremity->values[i]);
    }

    // leaf cell edge length = (V / num_leaf_cells) ^ 1/|P|
    index->settings->leaf_cell_edge_length = pow((index->settings->pivot_space_volume / index->settings->num_leaf_cells), (1.0/num_pivots));
    
    /* 
        Each leaf cell has a file called: 
        level_edge_center(mag,mean,stdev)_numVectors
        02_0.334_(0.33,12.282,7.008)
        level: leaf level id (m)
        edge length: edge length of the leaf cell
        mag: magnitude of the center vector of the leaf cell
        num vectors: number of vectors stored in leaf cell
        number of punctuation marks (underscores:3, parentheses:2, commas:2): total 7
    */
   
    float edge_length_size = ceil(log10(INT_MAX)+ 1);
    float center_vector_size = ceil(log10(INT_MAX)) * 3;
    float num_vectors_size = ceil(log10(SHRT_MAX) + 1);

    index->settings->max_filename_size = 1 + edge_length_size + center_vector_size +
                                            num_vectors_size + 7;

    // initialize file buffer manager
    if (!init_file_buffer_manager(index))
        exit_with_failure("Error in index.c: Could not initialize the file buffer manager for this index.");
    

    return OK;
}
/* get extremity vector of the pivot space, pivots must be outliers for extremity to be accurate */
vector * get_extremity(vector * pivot_vectors, unsigned int num_pivots)
{
    // extremity =  farthest vector in the pivot space (holds maximum distance to a pivot)
    vector * pivot_space_extremity = malloc(sizeof(struct vector));
    pivot_space_extremity->values = (v_type *) malloc(sizeof(v_type) * num_pivots);

    v_type max = FLT_MIN;
    for(int i = 0; i < num_pivots; i++)
    {
        for(int j = 0; j < num_pivots; j++)
            if(pivot_vectors[i].values[j] > max)
                max = pivot_vectors[i].values[j];
    }
    for(int i = 0; i < num_pivots; i++)
        pivot_space_extremity->values[i] = max;
        
    return pivot_space_extremity;
}

/* append vector to index */
response index_insert(pexeso_index * index, vector *vector)
{
    // get vector mapping in pivot space v -> v'
    struct vector * v_mapping = malloc(sizeof(struct vector));
    if(v_mapping == NULL)
        exit_with_failure("Error in index.c: Couldn't allocate memory for vector mapping.");
    v_mapping->values = malloc(sizeof(v_type) * index->settings->num_pivots);
    if(v_mapping->values == NULL)
        exit_with_failure("Error in index.c: Couldn't allocate memory for values of vector mapping.");
    
    map_vector(vector, index->settings->mtr_vector_length, v_mapping, index->settings->pivots_mtr, index->settings->num_pivots);

    printf("current vector, ");
    printf("in metric space :\n");
    print_vector(vector, index->settings->mtr_vector_length);
    printf("in pivot space :\n");
    print_vector(v_mapping, index->settings->num_pivots);

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
    

    printf("we reached a leaf cell with center :\n");
    print_vector(cell->center, index->settings->num_pivots);

    // append vector in metric format
    append_vector_to_cell(index, cell, vector);

    return OK;
}

void print_index(pexeso_index * index)
{
    printf("DISPLAY PEXESO INDEX...\n");
    printf("|       |       |       \n");
    printf("V       V       V       \n\n\n");
    level * level  = index->first_level;
    for(int i = 0; i < index->settings->num_levels; i++)
    {
        printf("Level %u:\n", level->id);
        printf("\tNumber of cells = %u\n", level->num_cells);
        printf("\tCell edge length = %f\n", level->cell_edge_length);
        printf("---------------------------------------------------------\n");
        for(int j = 0; j < level->num_cells; j++)
        {
            printf("++(Cell %d)\n--> Center vector: \t", j+1);
            print_vector(level->cells[j].center, index->settings->num_pivots);
            printf("\n--> Children: \n");
            if(level->cells[j].children == NULL)
                printf("(none)\n");
            else
            {
                for(int k = 0; k < level->cells[j].num_child_cells; k++)
                print_vector(level->cells[j].children[k].center, index->settings->num_pivots);
            }
        }
        printf("---------------------------------------------------------\n");
        
        level = level->next; 
    }
    printf("end of index.\n");
}