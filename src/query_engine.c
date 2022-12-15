#include <stdio.h>
#include <float.h>
#include "../include/hgrid.h"
#include "../include/cell.h"
#include "../include/match_map.h"
#include "../include/query_engine.h"
#include "../include/stats.h"
#include "../include/level.h"
/* verify candiate pairs */
enum response verify(struct grid *grid, struct pairs *pairs,
                     struct inv_index *index, struct match_map *match_map)
{
    float dist_threshold = grid->settings->query_settings->dist_threshold;
    float join_threshold = grid->settings->query_settings->join_threshold;
    unsigned int num_query_sets = grid->settings->query_settings->num_query_sets;
    unsigned int num_pivots = grid->settings->num_pivots;

    struct sid *query_set = malloc(sizeof(struct sid));
    query_set->table_id = 0;
    query_set->set_id = 0;
    query_set->set_size = 0;

    int map_idx = -1;

    // check if pairs list is not empty
    if (pairs->num_pairs == 0)
        exit_with_failure("Error in query_engine.c: Cannot verify 0 pairs! no matching pairs, no candidate pairs are found.");

    printf("\n\nProcess matching pairs...\n");
    if (pairs->matching_pairs != NULL)
    {
        // update match map for every set in a matching cell
        for (int i = 0; i < pairs->num_pairs; i++)
        {
            print_progress(i, pairs->num_pairs);

            struct matching_pair *mpair = &pairs->matching_pairs[i];
            struct vector *query_vector =  &pairs->query_vectors[i];
            
            query_set->table_id = pairs->query_vectors[i].table_id;
            query_set->set_id = pairs->query_vectors[i].set_id;
            query_set->set_size = pairs->query_vectors[i].set_size;
            
            map_idx = get_match_map_idx(match_map, num_query_sets, query_set);
            if(map_idx == -1)
                exit_with_failure("Error in match_map.c Couldn't find match map of query set.");
            
            RESET_QUERY_TIME()
            COUNT_QUERY_TIME_START
            if (pairs->has_matches[i])
            {
                for (int m = 0; m < mpair->num_match; m++) // loop through match cells
                {
                    struct cell *match_cell = mpair->cells[m];
                    COUNT_NEW_VISITED_CELL
                    COUNT_NEW_MATCH_CELL
                    if(match_cell->cell_size == 0)
                    {
                        printf("Warning in query_engine.c: empty matching cell!");
                        continue;
                    }
                    // get sets in matching cell, look for cell entry in inverted index
                    // int entry_idx = has_cell(index, match_cell, num_pivots);
                    // printf("\ncell, found index entry = %d, stored index entry = %d\n", entry_idx, match_cell->index_entry_pos);
                    int entry_idx = match_cell->index_entry_pos;
                    if (entry_idx == -1)
                    {
                        exit_with_failure("Error in query_engine.c: inverted index doesn't have entry for matching cell!");
                    }

                    struct entry *entry = &index->entries[entry_idx];
                    for (unsigned long s = 0; s < entry->num_sets; s++)
                    {
                        // update match count for all sets in matching cell
                        struct sid *curr_set = &index->distinct_sets[entry->sets[s]];
                        unsigned long set_idx = entry->sets[s];
                        if(match_map[map_idx].sets[set_idx].table_id != curr_set->table_id || match_map[map_idx].sets[set_idx].set_id != curr_set->set_id)
                            exit_with_failure("Error in query_engine.c: set position in match_map is diffrent from set position in inverted index.");
                    
                        update_match_count(match_map, map_idx, query_set, set_idx, join_threshold, query_vector->set_size);
                    }
                }
            }
            COUNT_QUERY_TIME_END
            COUNT_NEW_QUERY_TIME(query_time) // add curr query time to total query time (for all columns)
            match_map[map_idx].query_time += query_time;
        }
    }

    printf("\n\nProcess candidate pairs...\n");
    // <leaf cell -> vector tuples> array
    struct leaf_cell_vector_tuples * leaf_cells_tuples = init_leaf_cells_vectors_array(grid->settings);

    for (int i = 0; i < pairs->num_pairs; i++)
    {
        print_progress(i, pairs->num_pairs);

        struct candidate_pair *cpair = &pairs->candidate_pairs[i];
        struct vector *query_vector = &pairs->query_vectors[i]; // q'
        struct vector *query_vector_mtr = &pairs->query_vectors_mtr[i]; // q

        query_set->table_id = pairs->query_vectors[i].table_id;
        query_set->set_id = pairs->query_vectors[i].set_id;
        query_set->set_size = pairs->query_vectors[i].set_size;
        

        if (pairs->has_candidates[i]) // if current vector has candidate pair
        {
            map_idx = get_match_map_idx(match_map, num_query_sets, query_set);
            if(map_idx == -1)
                exit_with_failure("Error in match_map.c Couldn't find match map for query set.");
                
            RESET_QUERY_TIME()
            COUNT_QUERY_TIME_START

            for (int c = 0; c < cpair->num_candidates; c++) // loop through candidate cells
            {
                struct cell *candidate_cell = cpair->cells[c];
                COUNT_NEW_VISITED_CELL
                COUNT_NEW_CANDIDATE_CELL
                // get cell entry in the inverted index
                // int entry_idx = has_cell(index, candidate_cell, num_pivots);
                // printf("\ncell, found index entry = %d, stored index entry = %d\n", entry_idx, candidate_cell->index_entry_pos);
                int entry_idx = candidate_cell->index_entry_pos;
                if (entry_idx == -1)
                {
                    print_vector(candidate_cell->center, num_pivots);
                    exit_with_failure("Error in query_engine.c: inverted index doesn't have entry for candidate cell!");
                }
                
                //  get vectors stored in candidate cell
                struct vector_tuple *candidate_vectors = get_leaf_cell_vector_tuples(leaf_cells_tuples, candidate_cell,  grid->settings);
                if(candidate_vectors == NULL)
                    exit_with_failure("Error in query_engine.c: couldn't get leaf cell vector tuples.");

                // find entry of candidate cell in inverted index
                struct entry *entry = &index->entries[entry_idx];
                int curr_qvec_has_a_match = 0; // check if current query vector has  a match in current set or not
                unsigned long set_idx;

                // get vector (in metric and pivot space) in candidate cell
                for (int v = 0; v < candidate_cell->cell_size; v++)
                {
                    set_idx = candidate_vectors[v].mtr_vector->set_pos_in_inv_index;
                    // lemma 7: skip set if it cannot be joinable on "join_threshold" vectors
                    if ((query_set->set_size - match_map[map_idx].u[set_idx]) < ceil(join_threshold * query_set->set_size ))
                    {    
                        COUNT_USED_LEMMA(7)
                        continue;
                    }
                    int checked_vetors_in_ps = 0;
                    // printf("candidate vector: (%u, %u, %u)\n", candidate_vectors[v].mtr_vector->table_id, candidate_vectors[v].mtr_vector->set_id, candidate_vectors[v].mtr_vector->pos);
                    // if candidate vector can be filtered by lemma 1
                    if (pivot_filter(query_vector, candidate_vectors[v].ps_vector,
                                    grid->settings->num_pivots, dist_threshold))
                    {
                        COUNT_USED_LEMMA(1)
                        // printf("filterd by lemma 1\n");
                        // update mismatch count of curr set
                        update_mismatch_count(match_map, map_idx, set_idx);
                        checked_vetors_in_ps = 1;
                    }
                    // vector can be matched with q by lemma 2
                    else if (pivot_match(query_vector, candidate_vectors[v].ps_vector,
                                        grid->settings->num_pivots, dist_threshold))
                    {
                        COUNT_USED_LEMMA(2)
                        match_map[map_idx].has_match_for_curr_qvec[set_idx] = 1;
                        // printf("matched by lemma 2\n");
                        // update match count of curr set
                        update_match_count(match_map, map_idx, query_set, set_idx, join_threshold, query_vector->set_size);
                        checked_vetors_in_ps = 1;
                    }
                    else
                    {
                        COUNT_NEW_CHECKED_VECTOR_IN_MTR
                        // compute euclidean distance between v and q in metric space
                        float max_dist = 2.0, max_dist_no_sqrt = 4.0; // maximum distance between two normalized vectors
                        float dist_threshold_no_sqrt = (dist_threshold/max_dist) * max_dist_no_sqrt;
                        float d = euclidean_distance_cmp(candidate_vectors[v].mtr_vector,
                                                    query_vector_mtr, grid->settings->mtr_vector_length);
                        
                        // increase number of distance calculation
                        match_map[map_idx].num_dist_calc++;
                        
                        if (d <= dist_threshold_no_sqrt)
                        {
                            match_map[map_idx].has_match_for_curr_qvec[set_idx] = 1;
                            // printf("matched by euclidean dist\n");
                            update_match_count(match_map, map_idx, query_set, set_idx, join_threshold, query_vector->set_size);

                        }
                        else
                        {
                            // printf("filterd by euclidean dist\n");
                            update_mismatch_count(match_map, map_idx, set_idx);
                        }
                    }
                    if(checked_vetors_in_ps == 1)
                            COUNT_NEW_CHECKED_VECTOR_IN_PS

                    match_map[map_idx].total_checked_vectors++; 
                }   
            }
            // update |U| counter of vector in Q (column) that are not in S  (column) 
            update_zero_match_counter(&match_map[map_idx]);
            // for current set stop time and add it to query time
            COUNT_QUERY_TIME_END
            COUNT_NEW_QUERY_TIME(query_time) // add curr query time to total query time (for all columns)
            match_map[map_idx].query_time += query_time;
        }
    }
    // free memory
    free(query_set);
    if(!destroy_leaf_cells_vectors_array(leaf_cells_tuples, grid->settings))
        exit_with_failure("Error in query_engine.c: couldn't destroy leaf_cells_vectors_array!");
    return OK;
}

/* get candidate and matching cells of a query cell */
enum response block(struct cell *query_cell, struct cell *root_cell,
                    struct pairs *pairs, struct grid_settings *settings, struct match_map * match_map)
{
    float dist_threshold = settings->query_settings->dist_threshold;
    struct sid *query_set = malloc(sizeof(struct sid));
    query_set->table_id = 0;
    query_set->set_id = 0;
    query_set->set_size = 0;

    int map_idx = -1;
    
    is_empty(root_cell);
    is_empty(query_cell);
    printf("-b- block on cq %u, size = %ld\n", query_cell->id, query_cell->total_children_size);
    printf("-b- block on cr %u, size = %ld\n", root_cell->id, root_cell->total_children_size);
    for (int i = 0; i < query_cell->num_child_cells; i++)
    {
        struct cell *cq = &query_cell->children[i];
        if(is_empty(cq) == 1)
        {
            printf("\t-cq- cell %u (empty) ------------------------------------------------ \n", cq->id);
            continue;
        }
        else
            printf("\t-cq- cell %u (not empty) ------------------------------------------------ \n", cq->id);
        
        for (int j = 0; j < root_cell->num_child_cells; j++)
        {
            struct cell *cr = &root_cell->children[j];
            if(is_empty(cr) == 1)
                continue;
            printf("\t\t-cr- candidate cell %u\n", cr->id);
            // if cq and cr are leafs
            if (cq->is_leaf && cr->is_leaf)
            {
                if (cr->cell_size != 0) // cell are not empty
                {
                    // get all q' and q in cq
                    struct vector_tuple * query_vectors = get_vector_tuples(cq, settings, true);
                    // get vectors in cr cell
                    vector *cr_vectors = get_vectors_ps(cr, settings, false);

                    //for q, q' in cq
                    for (int q = 0; q < cq->cell_size; q++)
                    {
                        // check if current cell is already a candidate cell (added by quick browsing)
                        if(is_candidate_cell(pairs, query_vectors[q].ps_vector, cr, settings->num_pivots) == 1)
                        {
                            continue;
                        }
                        query_set->table_id = query_vectors[q].ps_vector->table_id;
                        query_set->set_id = query_vectors[q].ps_vector->set_id;
                        query_set->set_size = query_vectors[q].ps_vector->set_size;
                        
                        map_idx = get_match_map_idx(match_map, settings->query_settings->num_query_sets, query_set);
                        if(map_idx == -1)
                            exit_with_failure("Error in match_map.c Couldn't find match map of query set.");

                        RESET_QUERY_TIME()
                        COUNT_QUERY_TIME_START

                        // cr is a match of q' by lemma 5
                        if (vector_cell_match(query_vectors[q].ps_vector, cr, cr_vectors, settings, dist_threshold))
                        {
                            COUNT_USED_LEMMA(5)
                            if (!add_matching_pair(pairs, query_vectors[q].ps_vector, query_vectors[q].mtr_vector, cr, settings->num_pivots, settings->mtr_vector_length)) // add <q', cr>
                            {
                                exit_with_failure("Error in query_engine.c: couldn't add matching pair.");
                            }
                        }
                        else
                        {
                            // lemma 3
                            if (!vector_cell_filter(query_vectors[q].ps_vector, cr, cr_vectors, settings, dist_threshold))
                                add_candidate_pair(pairs, query_vectors[q].ps_vector, query_vectors[q].mtr_vector, cr, 
                                                    settings->num_pivots, settings->mtr_vector_length); // add <q', cr>
                            else
                            {
                                printf("-- fitered cell for vector %d -\n", q);
                                COUNT_NEW_FILTERED_CELL
                                COUNT_NEW_FILTERED_VECTORS(cr->cell_size)
                                COUNT_USED_LEMMA(3)
                            }
                        }
                        // for current set stop time and add it to query time
                        COUNT_QUERY_TIME_END
                        COUNT_NEW_QUERY_TIME(query_time) // add curr query time to total query time (for all columns)
                        match_map[map_idx].query_time += query_time;
                    }
                    // free memory
                    for(int q = 0; q < cq->cell_size; q++)
                    {
                        free(query_vectors[q].mtr_vector);
                        free(query_vectors[q].ps_vector);
                    }
                    free(query_vectors);

                    // free memory
                    for(int i = 0; i < cr->cell_size; i++)
                        free(cr_vectors[i].values);
                    free(cr_vectors);

                }
                else if (cq->cell_size == 0)
                {
                    // printf("Empty leaf query cell.\n\n");
                    break;
                }
            }
            else
            {
                // get vectors in query cell cq and vectors in cr cell
                // lemma 6
                if (cell_cell_match(cr, cq, settings, dist_threshold))
                {
                    COUNT_USED_LEMMA(6)
                    unsigned int num_cr_leaves = 0, max_leaf_idx = 0;
                    get_num_leaf_cells(cr, &num_cr_leaves);
                    max_leaf_idx = num_cr_leaves - 1;

                    struct cell **cr_leaves = NULL;
                    cr_leaves = malloc(sizeof(struct cell *) * (num_cr_leaves));
                    if (cr_leaves == NULL)
                        exit_with_failure("Error in cell.c: couldn't reallocate memory for root cell leaves.");

                    get_leaf_cells(cr, cr_leaves, &max_leaf_idx);

                    unsigned int num_cq_leaves = 0, max_qleaf_idx = 0;
                    get_num_leaf_cells(cq, &num_cq_leaves);
                    max_qleaf_idx = num_cq_leaves - 1;

                    struct cell **cq_leaves = NULL;
                    cq_leaves = malloc(sizeof(struct cell *) * (num_cq_leaves + 1));
                    if (cq_leaves == NULL)
                        exit_with_failure("Error in cell.c: couldn't reallocate memory for query cell leaves.");

                    get_leaf_cells(cq, cq_leaves, &max_qleaf_idx);

                    for (int ql = 0; ql < num_cq_leaves; ql++)
                    {
                        if(cq_leaves[ql]->cell_size == 0)
                            continue;

                        struct vector_tuple *query_vectors = get_vector_tuples(cq_leaves[ql], settings, true);
                        for (int q = 0; q < cq_leaves[ql]->cell_size; q++)
                        {
                            query_set->table_id = query_vectors[q].ps_vector->table_id;
                            query_set->set_id = query_vectors[q].ps_vector->set_id;
                            query_set->set_size = query_vectors[q].ps_vector->set_size;

                            map_idx = get_match_map_idx(match_map, settings->query_settings->num_query_sets, query_set);
                            if(map_idx == -1)
                                exit_with_failure("Error in match_map.c Couldn't find match map of query set.");

                            RESET_QUERY_TIME()
                            COUNT_QUERY_TIME_START

                            for (int l = 0; l < num_cr_leaves; l++)
                            {
                                // add cr to matching_pair
                                // printf("\n(+m) new match pair : (%.2f, %.2f)", cr_leaves[l]->center->values[0], cr_leaves[l]->center->values[1]);
                                if(cr_leaves[l]->cell_size != 0)
                                    add_matching_pair(pairs, query_vectors[q].ps_vector, query_vectors[q].mtr_vector, cr_leaves[l], settings->num_pivots, settings->mtr_vector_length);
                            }
                            // free memory
                            // free(query_vectors[q].mtr_vector->values);
                            free(query_vectors[q].mtr_vector);
                            // free(query_vectors[q].ps_vector->values);
                            free(query_vectors[q].ps_vector);

                            // for current set stop time and add it to query time
                            COUNT_QUERY_TIME_END
                            COUNT_NEW_QUERY_TIME(query_time) // add curr query time to total query time (for all columns)
                            match_map[map_idx].query_time += query_time;
                        }

                        // free memory
                        free(query_vectors);
                    }
                    // free memory
                    free(cr_leaves);
                    free(cq_leaves);
                }
                
                // lemma 4:
                else 
                {
                    if (cell_cell_filter(cr, cq, settings, dist_threshold))
                    {
                        printf("-f- fitered cell %u, size = %ld\n", cr->id, cr->total_children_size);
                        COUNT_NEW_FILTERED_CELL
                        COUNT_NEW_FILTERED_VECTORS(cr->total_children_size)
                        COUNT_USED_LEMMA(4)
                    }
                    else
                    {
                        // printf("-b- block on cell %u, size = %ld\n", cr->id, cr->total_children_size);
                        block(cq, cr, pairs, settings, match_map);
                    }
                }
            
            }
        }
    }
    free(query_set);
    return OK;
}

/* initialize query settings */
struct query_settings *init_query_settings(v_type dist_threshold, float join_threshold, int num_query_sets, int min_query_set_size, int max_query_set_size,
                                            double mtr_buffered_memory_size)
{
    struct query_settings *settings = malloc(sizeof(struct query_settings));
    settings->dist_threshold = dist_threshold;
    settings->join_threshold = join_threshold;
    settings->num_query_sets = num_query_sets;
    settings->min_query_set_size = min_query_set_size;
    settings->max_query_set_size = max_query_set_size;
    settings->mtr_buffered_memory_size = mtr_buffered_memory_size;
    return settings;
}
/* check if a vector is in the SQR of a query vector */
bool vector_in_SQR(struct vector * v, struct vector * q, unsigned int num_pivots, float dist_threshold)
{
    float max, min;
    for(int p = 0; p < num_pivots; p++)
    {
        max = q->values[p] + dist_threshold;
        min = q->values[p] - dist_threshold;
        if(v->values[p] <= max && v->values[p] >= min)
            continue;
        else
        {
            // printf("\n\nvector (%u, %u, %u) is NOT in SQR!\n\n",  v->table_id, v->set_id, v->pos);
            return 0;
        }
    }
    // printf("\n\nvector (%u, %u, %u) is in SQR!\n\n",  v->table_id, v->set_id, v->pos);
    return 1;
}

/* check if a vector is in the RGR of a query vector and pivot p */
bool vector_in_RQR(struct vector * v, struct vector * q, unsigned int p, float rqr)
{
    if(v->values[p] <= rqr)
    {
        // printf("\n\nvector (%u, %u, %u) is in RGR!\n\n",  v->table_id, v->set_id, v->pos);
        return 1;
    }
    return 0;
}

/*
    lemma 1:
    Given two vectors q and x, a set P of pivot vectors, a distance function d,
    and a threshold τ.
    if q matches x, then d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
    note: q and v are in pivot space.
    return OK if vector doesn't satisfy: d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
*/
enum response pivot_filter(struct vector *q, struct vector *x,
                           unsigned int num_pivots, float dist_threshold)
{
    for (int p = 0; p < num_pivots; p++)
    {
        // x match q means that: d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ i.e x cannot be filtered
        if (((q->values[p] - dist_threshold) <= x->values[p]) && ((q->values[p] + dist_threshold) >= x->values[p]))
            return FAILED;   
    }
    return OK;
}
/*
    lemma 2:
    Given two vectors q and x, a set P of pivot vectors, a distance function d,
    and a threshold τ, if there exists a pivot p ∈ P such that d(x, p) +d(q, p) ≤ τ,
    then q matches x.
    note: q and v are in pivot space.
*/
enum response pivot_match(struct vector *q, struct vector *x,
                          unsigned int num_pivots, v_type dist_threshold)
{
    // find pivot such that d(x, p) + d(q, p) ≤ τ
    for (int p = 0; p < num_pivots; p++)
    {
        if ((x->values[p] + q->values[p]) <= dist_threshold)
            return OK;
    }
    return FAILED;
}

/*
    lemma 3:
    Given a cell c and a mapped query vector q in the pivot space,
    if c ∩ SQR(q, τ) = ∅, then for any mapped vector x ∈ c,
    its original vector x does not match q.

    function returns OK if cell c can be filtered (pruned).
*/
enum response vector_cell_filter(struct vector *q, struct cell *cell, struct vector * cell_vectors,
                                 struct grid_settings * settings, float dist_threshold)
{
    // check if there is at least one vector in SQR of q
    for(int v = 0; v < cell->cell_size; v++)
    {
        if(vector_in_SQR(&cell_vectors[v], q, settings->num_pivots, dist_threshold))
        {
            for(int i = 0; i < cell->cell_size; i++)
                free(cell_vectors[i].values);
            free(cell_vectors);
            return FAILED;
        }
    }
    return OK;
}

/*
    lemma 4:
    Given a target cell c and a query cell cq in the pivot space,
    if c ∩ SQR(cq.center, (τ+cq.length / 2)) = ∅,
    then for any mapped vector x ∈ c and any query vector q ∈ cq,
    their original vectors do not match.

    function returns OK if cell c can be filtered (pruned).
*/
enum response cell_cell_filter(struct cell *cell, struct cell *query_cell,
                               struct grid_settings * settings, float dist_threshold)
{
    // check if there is at least one vector in SQR of cq.center
    float tau = dist_threshold + (query_cell->edge_length / 2);
    long unsigned int num_cell_vectors = 0;

    vector *cell_vectors = get_sub_cells_vectors_ps(cell, settings, &num_cell_vectors, false);

    for(int v = 0; v < num_cell_vectors; v++)
    {
        if(vector_in_SQR(&cell_vectors[v], query_cell->center, settings->num_pivots, tau))
        {
            // free memory
            for(int i = 0; i < num_cell_vectors; i++)
                free(cell_vectors[i].values);
            free(cell_vectors);

            return FAILED;
        }
    }

    // free memory
    for(int i = 0; i < num_cell_vectors; i++)
        free(cell_vectors[i].values);
    free(cell_vectors);
    return OK;
}


/*
    lemma 5:
    Given a target cell c and a mapped query vector q in the pivot space,
    if there exists a pivot p ∈ P such that c ∩ RQR(q, p, τ) = c,
    then for any vector x ∈ c, the original vector x matches the query vector q.

   function returns OK if cell c is a match to q.
*/
enum response vector_cell_match(struct vector *q, struct cell *cell, struct vector * cell_vectors,
                                struct grid_settings * settings, v_type dist_threshold)
{   
    float rqr;
    for (int p = 0; p < settings->num_pivots; p++)
    {
        rqr = dist_threshold - q->values[p];
        if(rqr < 0) // no rectangle query region for pivot p
            continue;

        bool  cell_in_rqr = 1; // assume cell is in rqr
        for(int v = 0; v < cell->cell_size; v++)
        {
            if(vector_in_RQR(&cell_vectors[v], q, p, rqr) == 0)
            {
                cell_in_rqr = 0; // cell is not in rqr (found one vector in c but not in rqr)
                break;  
            }
        }

        if(cell_in_rqr == 1)
        {
            // free memory
            for(int i = 0; i < cell->cell_size; i++)
                free(cell_vectors[i].values);
            free(cell_vectors);
            return OK;
        }
    }
    return FAILED;
}

/*
    lemma 6:
    Given a target cell c and a query cell cq in the pivot space,
    if there exists a pivot p ∈ P such that c ∩ min(RQR(q, p, τ)) = c,
    then for any mapped vector x ∈ c and any query vector q ∈ cq,
    their original vectors match.

    function returns OK if cell c is a match to cq.
*/
enum response cell_cell_match(struct cell *cell, struct cell *query_cell,
                              struct grid_settings * settings, float dist_threshold)
{
    long unsigned int num_cell_vectors = 0, num_query_vectors = 0;
    vector *cell_vectors, *query_vectors;

    if(cell->is_leaf == true && cell->cell_size != 0)
    {
        num_cell_vectors = cell->cell_size;
        cell_vectors = get_vectors_ps(cell, settings, false);
    }
    else if(cell->is_leaf == false)
        cell_vectors = get_sub_cells_vectors_ps(cell, settings, &num_cell_vectors, false);
    else
        return FAILED;

    if(query_cell->is_leaf == true && query_cell->cell_size != 0)
    {
        num_query_vectors = query_cell->cell_size;
        query_vectors = get_vectors_ps(query_cell, settings, true);
    }
    else if(query_cell->is_leaf == false)
        query_vectors = get_sub_cells_vectors_ps(query_cell, settings, &num_query_vectors, true);
    else
        return FAILED;
    

    float min_rqr;

    for(int q = 0; q < num_query_vectors; q++)
    {
        int p = min_RQR(&query_vectors[q], settings->num_pivots, dist_threshold);
       
        if(p == -1) // no rectangle query region for query vector q
            continue;

        min_rqr = dist_threshold - query_vectors[q].values[p];

        // if c ∩ min(RQR(q, p, τ)) = c
        bool  cell_in_rqr = 1; // assume cell is in rqr
        for(int v = 0; v < num_cell_vectors; v++)
        {
            if(vector_in_RQR(&cell_vectors[v], &query_vectors[q], p, min_rqr) == 0)
            {
                cell_in_rqr = 0; // cell is not in rqr (found one vector in c but not in rqr)
                break;  
            }
        }

        if(cell_in_rqr == 1)
        {
            // free memory
            for(int i = 0; i < num_cell_vectors; i++)
                free(cell_vectors[i].values);
            free(cell_vectors);

            for(int i = 0; i < num_query_vectors; i++)
                free(query_vectors[i].values);
            free(query_vectors);

            return OK;
        }
    }

    // free memory
    for(int i = 0; i < num_cell_vectors; i++)
        free(cell_vectors[i].values);
    free(cell_vectors);
    for(int i = 0; i < num_query_vectors; i++)
        free(query_vectors[i].values);
    free(query_vectors);

    return FAILED;
}

/* min rectagle query region RQR of a query in query_cell for pivot p, returns index of the pivot with the minimum rqr */
int min_RQR(struct vector *q, unsigned int num_pivots, float dist_threshold)
{
    float min_rqr = FLT_MAX, rqr;
    int bsf_p = -1;
    for (int p = 0; p < num_pivots; p++)
    {
        rqr = dist_threshold - q->values[p];
        if(rqr < 0)
            continue;
        // lead vectors stored in cell
        if (rqr < min_rqr)
        {
            min_rqr = rqr;
            bsf_p = p;
        }
    }

    if(min_rqr != FLT_MAX && bsf_p != -1)
        return bsf_p;

    return -1; // to indicate that no rqr is found
}
/* check if list of pairs already has query vector (todo: unvalid read of size 4)*/
int has_query_vector(struct pairs *pairs, vector *query_vector)
{
    if (pairs->num_pairs == 0)
        return -1;

    for (int i = 0; i < pairs->num_pairs; i++)
    {
        // printf("%p vs %p\n", pairs->query_vectors[i], query_vector);
        if (
            pairs->query_vectors[i].table_id == query_vector->table_id &&
            pairs->query_vectors[i].set_id == query_vector->set_id &&
            pairs->query_vectors[i].pos == query_vector->pos)
            return i;
    }

    return -1;
}

/* check if cell is already a candidate pair with query vector (through quick browsing) */
int is_candidate_cell(struct pairs * pairs, struct vector * q, struct cell * cell, unsigned int num_pivots)
{
    int q_idx = has_query_vector(pairs, q);

    // if query vector has no pairs
    if(q_idx == -1)
        return 0;

    struct candidate_pair *cp = &pairs->candidate_pairs[q_idx];

    struct cell * candidate_cell = NULL;
    int match = 1;

    for(int i = 0; i < cp->num_candidates; i++)
    {
        candidate_cell = cp->cells[i];
        match = 1;
        
        // for two cells to be the same they must have the same center vector
        if((cell->id == candidate_cell->id) && (cell->level->id == candidate_cell->level->id))
        {
            return 1;
        }
    }
    return 0;
}

/* add candidate pair */
enum response add_candidate_pair(struct pairs *pairs, struct vector *q, struct vector *q_mtr, struct cell *cell, unsigned int num_pivots, unsigned int mtr_vector_length)
{
    // printf("\n_____ add candidate! to (%u, %u, %u) \n", q->table_id, q->set_id, q->set_size);
    if (pairs == NULL)
        exit_with_failure("Error in query_engine.c: NULL pointer to result pairs!");

    unsigned int num_query_vectors = pairs->num_pairs;
    int q_idx = has_query_vector(pairs, q);

    // query vector doesn't exist in pairs list
    if (q_idx == -1 || num_query_vectors == 0)
    {
        // printf("\n\n\t*** *** *** Add query vector to pairs list (with its candidate)\n\n");
        if (q_idx == -1)
            q_idx = num_query_vectors; // add at last pos

        // alloc for new query vector
        pairs->query_vectors = realloc(pairs->query_vectors, sizeof(struct vector) * (num_query_vectors + 1));
        pairs->query_vectors_mtr = realloc(pairs->query_vectors_mtr, sizeof(struct vector) * (num_query_vectors + 1));
        pairs->candidate_pairs = realloc(pairs->candidate_pairs, sizeof(struct candidate_pair) * (num_query_vectors + 1));
        pairs->matching_pairs = realloc(pairs->matching_pairs, sizeof(struct matching_pair) * (num_query_vectors + 1));
        pairs->has_candidates = realloc(pairs->has_candidates, sizeof(bool) * (num_query_vectors + 1));
        pairs->has_matches = realloc(pairs->has_matches, sizeof(bool) * (num_query_vectors + 1));

        if
        (
            pairs->query_vectors == NULL || pairs->query_vectors_mtr == NULL || pairs->matching_pairs == NULL ||
            pairs->candidate_pairs == NULL || pairs->has_candidates == NULL || pairs->has_matches == NULL
        )            
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new match pair.");
        
        // copy vector
        pairs->query_vectors[q_idx].values = malloc(sizeof(v_type) * num_pivots);
        if (pairs->query_vectors[q_idx].values == NULL)
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new match pair pivot space values.");

        pairs->query_vectors_mtr[q_idx].values = malloc(sizeof(v_type) * mtr_vector_length);
        if (pairs->query_vectors_mtr[q_idx].values == NULL)
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new match pair metric space values.");


        for(int i = 0; i < num_pivots; i++)
            pairs->query_vectors[q_idx].values[i] = q->values[i];

        for(int i = 0; i < mtr_vector_length; i++)
            pairs->query_vectors_mtr[q_idx].values[i] = q_mtr->values[i];


        // table_id, set_ud and pos must be the same for q and q_mtr (same vector)
        pairs->query_vectors[q_idx].table_id = q->table_id;
        pairs->query_vectors[q_idx].set_id = q->set_id;
        pairs->query_vectors[q_idx].pos = q->pos;
        pairs->query_vectors[q_idx].set_size = q->set_size;

        pairs->query_vectors_mtr[q_idx].table_id = q_mtr->table_id;
        pairs->query_vectors_mtr[q_idx].set_id = q_mtr->set_id;
        pairs->query_vectors_mtr[q_idx].pos = q_mtr->pos;
        pairs->query_vectors_mtr[q_idx].set_size = q_mtr->set_size;

        pairs->candidate_pairs[q_idx].num_candidates = 0;
        pairs->candidate_pairs[q_idx].cells = NULL;
        pairs->has_candidates[q_idx] = true;

        pairs->has_matches[q_idx] = false;
        pairs->matching_pairs[q_idx].num_match = 0;
        pairs->matching_pairs[q_idx].cells = NULL;

        pairs->num_pairs++; // num distinct query vectors with pairs 
    }

    // add candidate 
    struct candidate_pair *cp = &pairs->candidate_pairs[q_idx];
    cp->cells = realloc(cp->cells, sizeof(struct cell *) * (cp->num_candidates + 1));
    if (cp->cells == NULL)
        exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new candidate pair cells.");
    cp->cells[cp->num_candidates] = cell;
    cp->num_candidates++;
    pairs->has_candidates[q_idx] = true;
    COUNT_ADD_CPAIR
    return OK;
}

/* add candidate pair */
enum response add_matching_pair(struct pairs *pairs, struct vector *q, struct vector *q_mtr, struct cell *cell, unsigned int num_pivots, unsigned int mtr_vector_length)
{
    if (pairs == NULL)
        exit_with_failure("Error in query_engine.c: NULL pointer to result pairs!");

    unsigned int num_query_vectors = pairs->num_pairs;
    int q_idx = has_query_vector(pairs, q);

    // query vector doesn't exist / the first to be added
    if (q_idx == -1 || num_query_vectors == 0)
    {
        // printf("\n\n\t*** *** *** Add query vector to pairs list (with its match)\n\n");
        if (q_idx == -1)
            q_idx = num_query_vectors; // add at last pos

        // alloc for new query vector
        pairs->query_vectors = realloc(pairs->query_vectors, sizeof(struct vector) * (num_query_vectors + 1));
        pairs->query_vectors_mtr = realloc(pairs->query_vectors_mtr, sizeof(struct vector) * (num_query_vectors + 1));
        pairs->candidate_pairs = realloc(pairs->candidate_pairs, sizeof(struct candidate_pair) * (num_query_vectors + 1));
        pairs->matching_pairs = realloc(pairs->matching_pairs, sizeof(struct matching_pair) * (num_query_vectors + 1));
        pairs->has_candidates = realloc(pairs->has_candidates, sizeof(bool) * (num_query_vectors + 1));
        pairs->has_matches = realloc(pairs->has_matches, sizeof(bool) * (num_query_vectors + 1));

        if
        (
            pairs->query_vectors == NULL || pairs->query_vectors_mtr == NULL || pairs->matching_pairs == NULL ||
            pairs->candidate_pairs == NULL || pairs->has_candidates == NULL || pairs->has_matches == NULL
        )            
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new match pair.");
        
        // copy vector
        pairs->query_vectors[q_idx].values = malloc(sizeof(v_type) * num_pivots);
        if (pairs->query_vectors[q_idx].values == NULL)
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new match pair pivot space values.");

        pairs->query_vectors_mtr[q_idx].values = malloc(sizeof(v_type) * mtr_vector_length);
        if (pairs->query_vectors_mtr[q_idx].values == NULL)
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new match pair metric space values.");


        for(int i = 0; i < num_pivots; i++)
            pairs->query_vectors[q_idx].values[i] = q->values[i];

        for(int i = 0; i < mtr_vector_length; i++)
            pairs->query_vectors_mtr[q_idx].values[i] = q_mtr->values[i];


        // table_id, set_ud and pos must be the same for q and q_mtr (same vector)
        pairs->query_vectors[q_idx].table_id = q->table_id;
        pairs->query_vectors[q_idx].set_id = q->set_id;
        pairs->query_vectors[q_idx].pos = q->pos;
        pairs->query_vectors[q_idx].set_size = q->set_size;

        pairs->query_vectors_mtr[q_idx].table_id = q_mtr->table_id;
        pairs->query_vectors_mtr[q_idx].set_id = q_mtr->set_id;
        pairs->query_vectors_mtr[q_idx].pos = q_mtr->pos;
        pairs->query_vectors_mtr[q_idx].set_size = q_mtr->set_size;
        
        pairs->matching_pairs[q_idx].num_match = 0;
        pairs->matching_pairs[q_idx].cells = NULL;
        pairs->has_matches[q_idx] = true;

        pairs->has_candidates[q_idx] = false;
        pairs->candidate_pairs[q_idx].num_candidates = 0;
        pairs->candidate_pairs[q_idx].cells = NULL;

        pairs->num_pairs++; // num distinct query vectors with pairs 
    }
    // add match
    struct matching_pair *mp = &pairs->matching_pairs[q_idx];
    mp->cells = realloc(mp->cells, sizeof(struct cell *) * (mp->num_match + 1));
    if (mp->cells == NULL)
        exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new candidate pair cells.");
    mp->cells[mp->num_match] = cell;
    mp->num_match++;
    pairs->has_matches[q_idx] = true;
    COUNT_ADD_MPAIR
    return OK;
}
/* destroy result pairs */
enum response destroy_pairs(struct pairs *pairs)
{
    if (pairs == NULL)
        exit_with_failure("Error in query_engine.c: NULL pointer to pairs!");

    if (pairs->num_pairs == 0)
    {
        free(pairs);
        return OK;
    }

    for (int i = pairs->num_pairs - 1; i >= 0; i--)
    {
        if (pairs->candidate_pairs != NULL)
        {
            if(pairs->has_candidates[i])
            {
                struct candidate_pair curr_cp = pairs->candidate_pairs[i];
                if(curr_cp.cells != NULL)
                    free(curr_cp.cells);
            }
        }

        if (pairs->matching_pairs != NULL)
        {
            
            if(pairs->has_matches[i])
            {
                struct matching_pair curr_mp = pairs->matching_pairs[i];
                if(curr_mp.cells != NULL)
                    free(curr_mp.cells);
                
            }
            
        }
        free(pairs->query_vectors[i].values);
        free(pairs->query_vectors_mtr[i].values);
    }
    if (pairs->matching_pairs != NULL)
        free(pairs->matching_pairs);
    if (pairs->candidate_pairs != NULL)
        free(pairs->candidate_pairs);
    if (pairs->query_vectors != NULL)
        free(pairs->query_vectors);
    if (pairs->query_vectors_mtr != NULL)
        free(pairs->query_vectors_mtr);
    if (pairs->has_candidates != NULL)
        free(pairs->has_candidates);
    if (pairs->has_matches != NULL)
        free(pairs->has_matches);
    free(pairs);
}

struct pairs *init_pairs()
{
    struct pairs *pairs = malloc(sizeof(struct pairs));
    if (pairs == NULL)
        exit_with_failure("Error in pexeso.c: Couldn't allocate memory for pairs!");

    pairs->num_pairs = 0u; // num vectors with candiate and/or matching pairs
    pairs->candidate_pairs = NULL;
    pairs->matching_pairs = NULL;
    pairs->query_vectors = NULL;
    pairs->query_vectors_mtr = NULL;
    pairs->has_candidates = false;
    pairs->has_matches = false;

    return pairs;
}