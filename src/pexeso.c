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
#include "../include/stats.h"
#include <string.h>
#include <valgrind/callgrind.h>

// run pexeso for one query column
enum response pexeso(struct sid * query_set, struct vector * query_vectors, unsigned int num_query_vectors,
                                struct grid * Dgrid, struct inv_index * inv_index)
{
    RESET_QUERY_TIME()
    COUNT_QUERY_TIME_START
    Dgrid->stats->total_queries_count++;
    reset_grid_stats(Dgrid);

    // build query grid and get list of query sets
    printf("\n\nBuilding query grid ... ");
    struct grid * Qgrid = (struct grid *) malloc(sizeof(struct grid));
    if (Qgrid == NULL)
        exit_with_failure("Error in pexeso.c: Couldn't allocate memory for grid!");

    // CALLGRIND_START_INSTRUMENTATION;
    // CALLGRIND_TOGGLE_COLLECT;

    if(!build_query_grid(Qgrid, Dgrid, inv_index, query_vectors, num_query_vectors))
        exit_with_failure("Error in pexeso.c: Couldn't build query grid!");

    printf("\n\nInit match maps... ");
    struct match_map * match_map =  init_match_map(inv_index, query_set);
    struct pairs * pairs = init_pairs(); // init list of candidate and matching pairs

    // quick browsing
    printf("\n\nQuick browsing... ");
    quick_browse(Dgrid, Qgrid, pairs, Dgrid->settings);
 
    // block (generate a list of candidate en matching pairs)
    printf("\n\nBlock...\n\n");
    block(Qgrid->root->cells, Dgrid->root->cells, pairs, Dgrid->settings, match_map);

    // verify (verify vectors in list of candidate pairs)
    printf("\n\nVerify...\n\n");
    verify(Dgrid, pairs, inv_index, match_map);

    // CALLGRIND_TOGGLE_COLLECT;
    // CALLGRIND_STOP_INSTRUMENTATION;

    COUNT_QUERY_TIME_END
    COUNT_NEW_QUERY_TIME(query_time) // add curr query time to total query time (for all columns)
    match_map->query_time += query_time;

    // print Query grid
    // dump_grid_to_console(Qgrid);

    // print match map after quering
    dump_csv_results_to_console(match_map, 0);
    
    if(!save_results_to_disk(Dgrid, Qgrid, match_map))
        exit_with_failure("Error in pexeso.c: Couldn't save query results to disk.");


    Dgrid->stats->total_query_time += total_query_time;
    Dgrid->stats->loaded_query_files_count += loaded_query_files_count;
    Dgrid->stats->loaded_query_sets_count += loaded_query_sets_count;
    Dgrid->stats->loaded_query_files_size += loaded_query_files_size;
    Dgrid->stats->loaded_qvec_count += loaded_qvec_count;
    for(int i = 0; i < 7; i++)
        Dgrid->stats->used_lemmas_count[i] = used_lemmas_count[i];

    Dgrid->stats->filtered_cells_count = filtered_cells_count;
    Dgrid->stats->visited_cells_count = visited_cells_count;
    Dgrid->stats->visited_matching_cells_count = visited_matching_cells_count;
    Dgrid->stats->visited_candidate_cells_count = visited_candidate_cells_count;

    Dgrid->stats->filterd_vectors_count = filterd_vectors_count;
    Dgrid->stats->checked_vectors_in_ps_count = checked_vectors_in_ps_count;
    Dgrid->stats->checked_vectors_in_mtr_count = checked_vectors_in_mtr_count;

    Dgrid->stats->count_add_cpair = count_add_cpair;
    Dgrid->stats->count_add_mpair = count_add_mpair;
    Dgrid->stats->count_dist_calc = count_dist_calc;

    /* end of endexing and quering */
    COUNT_TOTAL_TIME_END
    Dgrid->stats->total_time = total_time;

    /* print grid statistics */
    print_grid_stats(Dgrid);

    /* destroy match map */
    if(!match_maps_destroy(match_map, 1))
        exit_with_failure("Error main.c: Couldn't destroy match map.\n");
    
    // free memory: destroy query grid and result pairs
    destroy_pairs(pairs);
    query_grid_destroy(Qgrid);
}

/* quick browsing: evaluate leaf cells in Qgrid in Dgrid inverted index and get candidate pairs */
void quick_browse(struct grid * Dgrid, struct grid * Qgrid, struct pairs * pairs, struct grid_settings * settings)
{
    unsigned int num_pivots = Dgrid->settings->num_pivots;

    struct level * level = Dgrid->root;
    struct level * qlevel = Qgrid->root; 

    while(!level->is_leaf) // go to leaf level in data grid
        level = level->next;

    while(!qlevel->is_leaf) // go to leaf level in query grid
        qlevel = qlevel->next;

    struct cell * cr, *cq;
    for(int i = 0; i < level->num_cells; i++)
    {
        cr = &level->cells[i];
        if(cr->cell_size == 0) continue; // skip if cr cell is empty
        for(int j = 0; j < qlevel->num_cells; j++)
        {
            cq = &qlevel->cells[j];
            if(cq->cell_size == 0) continue; // skip if query cell is empty

            if(cr->id == cq->id) // cq and cr refer to the same region
            {
                // create candiate pair <q', cr> for every query vector in cq
                // get all q' and q in cq
                struct vector_tuple * query_vectors = get_vector_tuples(cq, settings, true);
                for (int q = 0; q < cq->cell_size; q++)
                {
                    // add <q', cr>
                    add_candidate_pair(pairs, query_vectors[q].ps_vector, query_vectors[q].mtr_vector, cr, settings->num_pivots, settings->mtr_vector_length); 
                    
                    // free(query_vectors[q].mtr_vector->values);
                    free(query_vectors[q].mtr_vector);
                    // free(query_vectors[q].ps_vector->values);
                    free(query_vectors[q].ps_vector);
                }
                // free memory
                free(query_vectors);
            }
            
        }
    }
}

/* create one query grid  */
enum response build_query_grid(struct grid * Qgrid, struct grid * Dgrid, inv_index * inv_index, struct vector * query_vectors, unsigned int num_query_vectors)
{   
    long unsigned int num_query_files = 0; // always = 1;
    unsigned int mtr_query_vector_length = Dgrid->settings->mtr_vector_length;
    unsigned int base = Dgrid->settings->base;

    unsigned int mtr_buffered_memory_size = Dgrid->settings->query_settings->mtr_buffered_memory_size;
    unsigned int ps_buffered_memory_size = Dgrid->settings->query_settings->ps_buffered_memory_size;
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

    // printf("\n\n\n\t..................................\n");
    // printf("\t::      BUILD QUERY INDEX       ::\n");
    // printf("\t..................................\n\n\n");


    /* pivot space extremity */
    // printf("\n\n(!) Use pivot vectors in Dgrid ...Extremity vector (in pivot space):\n");
    vector * pivot_space_extremity = Dgrid->settings->pivot_space_extremity;
    vector * pivots_mtr = Dgrid->settings->pivots_mtr; // pivot vectors (in metric space)
    vector * pivots_ps = Dgrid->settings->pivots_ps; // pivot vectors (in pivot space)
    // print_vector(pivot_space_extremity, num_pivots);

    /* initialize grid */
    if (!init_grid(Dgrid->settings->work_directory, 0, num_pivots, pivots_mtr, pivots_ps, pivot_space_extremity, 
                    num_levels, num_query_vectors, base, mtr_query_vector_length, 
                    mtr_buffered_memory_size, max_leaf_size, track_vector, 
                    true, Dgrid->settings->query_settings, Qgrid, 0))
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
    
    for(int qv = 0; qv < num_query_vectors; qv++)
    {
        if(!grid_insert(Qgrid, inv_index, &query_vectors[qv]))
            exit_with_failure("Error in pexeso.c: Couldn't insert vector in query grid.");
    }

    return OK;
}