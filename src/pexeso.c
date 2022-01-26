#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/gsl_matrix.h"
#include "../include/select_pivots.h"
#include "../include/match_map.h"
#include "../include/query_engine.h"
#include "../include/file_loader.h"
#include "../include/pexeso.h"

/* pexeso set similarity search algorithm (result in match map of Dgrid) */
void pexeso(const char * query_file_dir, struct grid * Dgrid, struct inv_index * inv_index,
             struct match_map * match_map)
{
    long long unsigned int num_query_vectors = 0;
    // make grid for query file
    struct grid * Qgrid = make_query_grid(Dgrid, inv_index, query_file_dir, &num_query_vectors);

    // (todo) quick browsing 

    // init list of candidate and matching pairs
    struct pairs * pairs = init_pairs();
    
    //block
    block(Qgrid->root->cells, Dgrid->root->cells, pairs, Dgrid->settings);

    //verify
    // verify(Dgrid, pairs, inv_index, match_map, Qgrid->total_records);

    // print Query grid
    // dump_grid_to_console(Qgrid);

    // print match map after quering
    dump_match_map_to_console(match_map);

    // free memory: destroy query grid and result pairs
    destroy_pairs(pairs);
    query_grid_destroy(Qgrid);
}

/* quick browsing: evaluate leaf cells in Qgrid in Dgrid inverted index and get candidate pairs*/
void quick_browse(struct grid * Dgrid, struct grid * Qgrid)
{
    struct cell * cell = Qgrid->root->cells; // get cell in root level
    struct cell ** qgrid_leaves = NULL;
    unsigned int  * num_leaves  = 0;
    get_leaf_cells(cell, qgrid_leaves, num_leaves);
    
    for(int l = 0; l < *num_leaves; l++)
    {
        
    }
}
/* create query grid  */
struct grid * make_query_grid(struct grid * Dgrid, inv_index * inv_index, const char * query_file_dir, long long unsigned int * num_query_vectors)
{   
    long unsigned int num_query_files = 0; // always = 1;
    unsigned int mtr_query_vector_length = Dgrid->settings->mtr_vector_length;
    unsigned int base = Dgrid->settings->base;

    unsigned int mtr_buffered_memory_size = 1;
    unsigned int ps_buffered_memory_size = 1;
    unsigned int max_leaf_size = Dgrid->settings->max_leaf_size; // max vectors in one leaf cell
    char * query_grid_dir = "/home/jaouhara/Projects/pexeso-debbug/pexeso/Qgrid/";

    unsigned int track_vector = 1; // track vectors id (table_id, column_id)
    
    /* grid settings */
    unsigned int num_levels = Dgrid->settings->num_levels;  // Dgrid and Qgrid are constructed with the same level number
    unsigned int num_pivots = Dgrid->settings->num_pivots;  // number of pivots
    unsigned int fft_scale = 13;   // constant for finding |P| * fft_scale candidate pivots, a good choice of fft_scale is approximately 30 (in paper) and 13 with experiments.

    /* query settings */
    unsigned int join_threshold = Dgrid->settings->query_settings->join_threshold; // T 
    v_type dist_threshold = Dgrid->settings->query_settings->dist_threshold; // tau

    printf("\n\n\n\t..................................\n");
    printf("\t::      BUILD QUERY INDEX       ::\n");
    printf("\t..................................\n\n\n");

    /* read all vectors in the data set */
    printf("\n\n\n\nReading query info...");
    get_dataset_info(query_file_dir, &num_query_files, num_query_vectors, &mtr_query_vector_length);
    printf("(OK)\n");
    printf("\tNumber of (.bin) files = %lu\n\tNumber of vectors = %llu\n\t"
            "Vector length in mtric spaces = %u\n\n", num_query_files, *num_query_vectors,
             mtr_query_vector_length);

    // worst case: senario all query vectors will end up in one leaf cell
    max_leaf_size = *num_query_vectors;


    /* pivot space extremity */
    printf("\n\n(!) Use pivot vectors in Dgrid ...Extremity vector (in pivot space):\n");
    
    // vector * pivot_space_extremity = get_extremity(pivots_ps, num_pivots);
    vector * pivot_space_extremity = Dgrid->settings->pivot_space_extremity;
    vector * pivots_mtr = Dgrid->settings->pivots_mtr; // pivot vectors (in metric space)
    vector * pivots_ps = Dgrid->settings->pivots_ps; // pivot vectors (in pivot space)
    print_vector(pivot_space_extremity, num_pivots);
    
    printf("(OK)\n");

    /* initialize grid */
    printf("\n\nInitialize grid... ");
    struct grid * Qgrid = (struct grid *) malloc(sizeof(struct grid));
    if (Qgrid == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for grid!");

    if (!init_grid(query_grid_dir, num_pivots, pivots_mtr, pivots_ps, pivot_space_extremity, 
                    num_levels, *num_query_vectors, base, mtr_query_vector_length, 
                    mtr_buffered_memory_size, ps_buffered_memory_size, max_leaf_size, track_vector, 
                    true, Dgrid->settings->query_settings, Qgrid))
        exit_with_failure("Error in main.c: Couldn't initialize grid!");

    printf("(OK)\n");


    /* Build levels */
    printf("\n\nBuild query grid levels... ");
    if(!init_root(Qgrid))
        exit_with_failure("Error in main.c: Couldn't initialize root level!");

    if (!init_first_level(Qgrid))
        exit_with_failure("Error in main.c: Couldn't initialize first level!");

    if (!init_levels(Qgrid))
        exit_with_failure("Error in main.c: Couldn't initialize grid levels!");
    printf("(OK)\n");


    /* insert dataset in grid (read and index all vectors in the data set) */
    printf("Index dataset vectors and build inverted index...");
    struct inv_index * index = malloc(sizeof(struct inv_index));
    index->num_entries = 0;
    if(index == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for inverted index.");

    if (!index_binary_files(Qgrid, inv_index, query_file_dir, num_query_files, base))
        exit_with_failure("Error in main.c: Something went wrong, couldn't index binary files.");
    printf("(OK)\n");

    /* print grid */
    // dump_grid_to_console(Qgrid);

    // (todo) index query vectors without creating an inverted index
    /* destroy inverted index (inverted index not required for queries*/
    if(!inv_index_destroy(index))
        exit_with_failure("Error main.c: Couldn't destroy inverted index.\n");

    return Qgrid;
}