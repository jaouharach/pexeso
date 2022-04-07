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
#include <string.h>

/* pexeso set similarity search algorithm (result in match map of Dgrid) */
void pexeso(const char * query_file_dir, struct grid * Dgrid, struct inv_index * inv_index)
{
    long long unsigned int num_query_vectors = 0;
    unsigned int num_query_sets = 0;

    // build query grid and get list of query sets
    printf("\n\nBuilding query grid... ");
    struct grid * Qgrid = (struct grid *) malloc(sizeof(struct grid));
    if (Qgrid == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for grid!");

    struct sid * query_sets = build_query_grid(Qgrid, Dgrid, inv_index, query_file_dir, &num_query_vectors);

    num_query_sets = Dgrid->settings->query_settings->num_query_sets;

    //initialize match maps for query sets
    printf("\n\nMaking match maps... ");
    struct match_map * match_map =  init_match_maps(inv_index, query_sets, num_query_sets);


    // (todo) quick browsing 

    // init list of candidate and matching pairs
    struct pairs * pairs = init_pairs();
    
    // block (generate a list of candidate en matching pairs)
    block(Qgrid->root->cells, Dgrid->root->cells, pairs, Dgrid->settings, match_map);

    // verify (verify vectors in list of candidate pairs)
    verify(Dgrid, pairs, inv_index, match_map);

    // print Query grid
    // dump_grid_to_console(Qgrid);

    // print match map after quering
    for(int m = 0; m < num_query_sets; m++)
        dump_match_map_to_console(match_map, m);
    
    if(!save_results_to_disk(Dgrid, Qgrid, match_map))
        exit_with_failure("Error in pexeso.c: Couldn't save query results to disk.");


    /* destroy match map */
    if(!match_maps_destroy(match_map, num_query_sets))
        exit_with_failure("Error main.c: Couldn't destroy match map.\n");
    
    // free memory: destroy query grid and result pairs
    destroy_pairs(pairs);
    query_grid_destroy(Qgrid);
}

/* quick browsing: evaluate leaf cells in Qgrid in Dgrid inverted index and get candidate pairs */
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
struct sid * build_query_grid(struct grid * Qgrid, struct grid * Dgrid, inv_index * inv_index, const char * query_file_dir, long long unsigned int * num_query_vectors)
{   
    long unsigned int num_query_files = 0; // always = 1;
    unsigned int mtr_query_vector_length = Dgrid->settings->mtr_vector_length;
    unsigned int base = Dgrid->settings->base;

    unsigned int mtr_buffered_memory_size = 1;
    unsigned int ps_buffered_memory_size = 1;
    unsigned int max_leaf_size = Dgrid->settings->max_leaf_size; // max vectors in one leaf cell

    unsigned int track_vector = 1; // track vectors id (table_id, column_id)
    
    /* grid settings */
    unsigned int num_levels = Dgrid->settings->num_levels;  // Dgrid and Qgrid are constructed with the same level number
    unsigned int num_pivots = Dgrid->settings->num_pivots;  // number of pivots
    unsigned int fft_scale = 13;   // constant for finding |P| * fft_scale candidate pivots, a good choice of fft_scale is approximately 30 (in paper) and 13 with experiments.

    /* query settings */
    float join_threshold = Dgrid->settings->query_settings->join_threshold; // T 
    float dist_threshold = Dgrid->settings->query_settings->dist_threshold; // tau
    int num_query_sets = Dgrid->settings->query_settings->num_query_sets;
    int min_query_set_size = Dgrid->settings->query_settings->min_query_set_size;
    int max_query_set_size = Dgrid->settings->query_settings->max_query_set_size;

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
    if (!init_grid(Dgrid->settings->work_directory, num_pivots, pivots_mtr, pivots_ps, pivot_space_extremity, 
                    num_levels, *num_query_vectors, base, mtr_query_vector_length, 
                    mtr_buffered_memory_size, ps_buffered_memory_size, max_leaf_size, track_vector, 
                    true, Dgrid->settings->query_settings, Qgrid))
        exit_with_failure("Error in main.c: Couldn't initialize grid!");

    // to allow getting query vectors from disk  using Dgrid settings
    Dgrid->settings->query_root_directory = Qgrid->settings->root_directory;

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
    printf("Index query vectors...");
    
    
    struct sid * query_sets = index_query_binary_files(Qgrid, Dgrid, inv_index, query_file_dir, num_query_files, base, min_query_set_size, max_query_set_size);
    if (query_sets == NULL)
        exit_with_failure("Error in main.c: Something went wrong, couldn't get list of query sets after indexing binary files.");
    printf("(OK)\n");

    return query_sets;
}