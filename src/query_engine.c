#include <stdio.h>
#include <float.h>
#include "../include/hgrid.h"
#include "../include/cell.h"
#include "../include/match_map.h"
#include "../include/query_engine.h"

// /* verify candiate pairs */
// enum response verify(struct grid * grid, struct matching_pair * mpair, struct candidate_pair * cpair,
//             struct inv_index * index, struct match_map * match_map, unsigned int query_set_size)
// {   
//     v_type dist_threshold = grid->settings->query_settings->dist_threshold;
//     unsigned int join_threshold = grid->settings->query_settings->join_threshold;

//     // update match map for every set in a matching cell
//     for(int m = 0; m < mpair->num_match;  m++)
//     {
//         struct cell * match_cell = mpair->cells[m];

//         // get entry of matching cell in inverted index
//         int entry_idx = has_cell(index, match_cell);
//         if(entry_idx == -1)
//             exit_with_failure("Error in query_engine.c: inverted index doesn't has entry fro matching cell!");

//         struct entry * curr_entry = &index->entries[entry_idx];
//         for(int s = 0; s < curr_entry->num_sets; s++)
//         {
//             //update match count for all sets in matching cell
//             struct sid * curr_set = &index->distinct_sets[curr_entry->sets[s]];
//             printf("\nmatch ++ !!!\n");
//             update_match_count(match_map, curr_set);
//         }
//     }

//     // verify candidates
//     for(int c = 0; c < cpair->num_candidates;  c++)
//     {
//         struct cell * candidate_cell = cpair->cells[c];
//         struct vector * query_vector = cpair->query_vectors[c]; // in ps
        
//         int entry_idx = has_cell(index, candidate_cell);
//         if(entry_idx == -1)
//             exit_with_failure("Error in query_engine.c: inverted index doesn't has entry fro matching cell!");

//         // find entry of candidate cell in inverted index
//         struct entry * candidate_cell_entry = &index->entries[entry_idx];

//         // get sets (columns) in candidate_cell
//         // for every set in candidate cell
//         for(int s = 0; s < candidate_cell_entry->num_sets; s++)
//         {
//             struct sid * curr_set = &index->distinct_sets[candidate_cell_entry->sets[s]];

//             if(curr_set == NULL)
//                 exit_with_failure("Error in query_engine.c: NULL set in candidate entry!");
            
//             // set id position in match map
//             int set_idx = has_set(match_map, curr_set);
//             // lemma 7: skip set if it cannot be joinable on "join_threshold" vectors
//             if(query_set_size - match_map->mismatch_count[set_idx] < join_threshold)
//                 continue;

//             // get vector (in metric and pivot space) in candidate cell (in metric space)
//             struct vector_tuple * candidate_vectors = get_vector_tuples(candidate_cell, grid->settings);
//             for(int v = 0; v < candidate_cell->cell_size; v++)
//             {
//                 // if vector belongs to set (curr set)
//                 if (candidate_vectors[v].mtr_vector->set_id == curr_set->set_pos 
//                 && candidate_vectors[v].mtr_vector->table_id == curr_set->table_id)
//                 {
//                     // if candidate vector can be filtered by lemma 1
//                     if(pivot_filter(query_vector, candidate_vectors[v].ps_vector, 
//                     grid->settings->num_pivots, dist_threshold))
//                     {
//                         // update mismatch count of curr set
//                         printf("\nmiss match ++ !!!\n");
//                         match_map->mismatch_count[set_idx]++;
//                     }

//                     // vector can be matched with q by lemma 2
//                     else if (pivot_match(query_vector, candidate_vectors[v].ps_vector, 
//                     grid->settings->num_pivots, dist_threshold))
//                     {
//                         // update match count of curr set
//                         printf("\nmatch ++ !!!\n");
//                         match_map->match_count[set_idx]++;
//                     }
//                     else
//                     {
//                         // compute euclidean distance between v and q in metric space
//                         v_type d = euclidean_distance(candidate_vectors[v].mtr_vector, 
//                         candidate_vectors[v].mtr_vector, grid->settings->mtr_vector_length);
//                         if(d <= dist_threshold)
//                         {
//                             printf("\nmatch ++ !!!\n");
//                             match_map->match_count[set_idx]++;
//                         }   
//                         else
//                         {
//                             printf("\nmiss match ++ !!!\n");
//                             match_map->mismatch_count[set_idx]++;
//                         }
//                     }
//                 }
//                 if(match_map->match_count[set_idx] > join_threshold)
//                 {
//                     // mark current set (curr_entry) as joinable
//                     match_map->joinable[set_idx] = true;
//                     continue;
//                 }
//             }
//         }
//     }

//     return OK;
// }

/* get candidate and matching cells of a query cell */
enum response block(struct cell *query_cell, struct cell * root_cell, 
                    struct pairs * pairs, struct grid_settings * settings)
{
    v_type dist_threshold = settings->query_settings->dist_threshold;
    //forech child cell of query cell
    for(int i = 0; i < query_cell->num_child_cells; i++)
    {
        struct cell * cq = &query_cell->children[i];
        //forech child cell of r cell
        for(int j = 0; j < root_cell->num_child_cells; j++)
        {
            struct cell * cr = &root_cell->children[j];

            //if cq and cr are leafs
            if(cq->is_leaf && cr->is_leaf)
            {
                if(cq->cell_size != 0 && cr->cell_size != 0) // cell are not empty
                {   
                    // query vector in pivot space
                    struct vector ** query_vector = get_vectors_ps(cq, settings->mtr_vector_length);
                    
                    for(int q = 0; q < cq->cell_size; q++)
                    {   
                        printf("(%u, %u, %u)\n", query_vector[q]->table_id, query_vector[q]->set_id, query_vector[q]->pos);
                        // cr is a match of q by lemma 5
                        if(vector_cell_match(query_vector[q], cr, settings->num_pivots, dist_threshold))
                        {
                            // add cr to matching_pair
                            add_matching_pair(pairs, query_vector[q], cr);
                        }
                        // if cr cannot be pruned by lemma 3
                        else if(!vector_cell_filter(query_vector[q], cr, settings->num_pivots, dist_threshold))
                        {
                            // add cr to candidate_pair
                            add_candidate_pair(pairs, query_vector[q], cr);
                        }
                    }
                    // free memory
                    free(query_vector);
                }
            }
            else
            {
                // lemma 6
                if(cell_cell_match(cq, cr, settings->num_pivots, dist_threshold))
                {   
                    unsigned int num_cr_leaves  = 0;
                    get_num_leaf_cells(cr, &num_cr_leaves);

                    struct cell ** cr_leaves = NULL;
                    cr_leaves = malloc(sizeof(struct cell *) * (num_cr_leaves + 1));
                    if(cr_leaves == NULL)
                        exit_with_failure("Error in cell.c: couldn't reallocate memory for root cell leaves.");
                    
                    get_leaf_cells(cr, cr_leaves, &num_cr_leaves);

                    
                    unsigned int num_cq_leaves  = 0;
                    get_num_leaf_cells(cq, &num_cq_leaves);

                    struct cell ** cq_leaves = NULL;
                    cq_leaves = malloc(sizeof(struct cell *) * (num_cq_leaves + 1));
                    if(cq_leaves == NULL)
                        exit_with_failure("Error in cell.c: couldn't reallocate memory for query cell leaves.");
                    
                    get_leaf_cells(cq, cq_leaves, &num_cq_leaves);

                    for(int ql = 0; ql < num_cq_leaves; ql++)
                    {
                        struct vector ** query_vector = get_vectors_ps(cq_leaves[ql], settings->num_pivots);
                        for(int q = 0; q < cq->cell_size; q++)
                        {   
                            for(int l = 0; l < num_cr_leaves; l++)
                                // add cr to matching_pair
                                add_matching_pair(pairs, query_vector[q], cr_leaves[l]);
                                
                        }
                        // free memory
                        free(query_vector);
                    }
                }
                // lemma 4: cr cannot be pruned
                else if (!cell_cell_filter(cq, cr, settings->num_pivots, dist_threshold))
                {
                    block(cq, cr, pairs, settings);
                }
            }
        }
    }
    return OK;
}

/* initialize query settings */
struct query_settings * init_query_settings(v_type dist_threshold, unsigned int join_threshold)
{
    struct query_settings * settings = malloc(sizeof(struct query_settings));
    settings->dist_threshold = dist_threshold;
    settings->join_threshold = join_threshold;

    return settings;
}
/* 
    Given two vectors q and x, a set P of pivot vectors, a distance function d, 
    and a threshold τ.
    if q matches x, then d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
    note: q and v are in pivot space.
    return OK if vector doesn't satisfy: d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
*/
enum response pivot_filter(struct vector * q, struct vector * x, 
                            unsigned int num_pivots, v_type dist_threshold)
{
    // find pivot such that d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
    bool filter = true;
    for(int p = 0; p < num_pivots; p++)
    {
        if(
            ((q->values[p] - dist_threshold) <= x->values[p])
            &&
            ((q->values[p] + dist_threshold) >= x->values[p])
        )
            filter = false;
    }
    return filter ? OK : FAILED;
}
/* 
    Given two vectors q and x, a set P of pivot vectors, a distance function d, 
    and a threshold τ, if there exists a pivot p ∈ P such that d(x, p) +d(q, p) ≤ τ,
    then q matches x.
    note: q and v are in pivot space.
*/
enum response pivot_match(struct vector * q, struct vector * x, 
                            unsigned int num_pivots, v_type dist_threshold)
{
    // find pivot such that d(x, p) +d(q, p) ≤ τ
    for(int p = 0; p < num_pivots; p++)
    {
        if((x->values[p] + q->values[p]) <= dist_threshold)
            return OK;
    }
    return FAILED;
}

/* 
   Given a cell c and a mapped query vector q in the pivot space,
   if c ∩ SQR(q, τ) = ∅, then for any mapped vector x ∈ c, 
   its original vector x does not match q.

   function returns OK if cell c can be filtered (pruned).
*/
enum response vector_cell_filter(struct vector * q, struct cell * cell, 
                            unsigned int num_pivots, v_type dist_threshold)
{
    if(euclidean_distance(cell->center, q, num_pivots) >= dist_threshold)
        return OK;

    return FAILED;
}

/* 
   Given a target cell c and a query cell cq in the pivot space, 
   if c ∩ SQR(cq.center, (τ+cq.length / 2)) = ∅, 
   then for any mapped vector x ∈ c and any query vector q ∈ cq, 
   their original vectors do not match.

   function returns OK if cell c can be filtered (pruned).
*/
enum response cell_cell_filter(struct cell * cell, struct cell * query_cell, 
                            unsigned int num_pivots, v_type dist_threshold)
{
    if(euclidean_distance(cell->center, query_cell->center, num_pivots)
                            >= ((dist_threshold + query_cell->edge_length) / 2))
        return OK;
        
    return FAILED;
}

/* 
    Given a target cell c and a mapped query vector q in the pivot space,
    if there exists a pivot p ∈ P such that c ∩ RQR(q, p, τ) = c, 
    then for any vector x ∈ c, the original vector x matches the query vector q.

   function returns OK if cell c is a match to q.
*/
enum response vector_cell_match(struct vector * q, struct cell * cell, 
                            unsigned int num_pivots, v_type dist_threshold)
{
    for(int p = 0; p < num_pivots; p++)
    {
        if((cell->center->values[p] * 2) - dist_threshold + q->values[p] == 0)
            return OK;
    }
    return FAILED;
}

/* 
   Given a target cell c and a query cell cq in the pivot space, 
   if there exists a pivot p ∈ P such that c ∩ min(RQR(q, p, τ)) = c, 
   then for any mapped vector x ∈ c and any query vector q ∈ cq,
   their original vectors match.

   function returns OK if cell c is a match to cq.
*/
enum response cell_cell_match(struct cell * cell, struct cell * query_cell, 
                            unsigned int num_pivots, v_type dist_threshold)
{
    // min RQR
    v_type min_rqr;

    for(int p = 0; p < num_pivots; p++)
    {
        min_rqr = min_RQR(query_cell, num_pivots, p, dist_threshold);
        // if c ∩ min(RQR(q, p, τ)) = c
        if((cell->center->values[p] * 2) - min_rqr == 0) 
            return OK;
    }
        
    return FAILED;
}

/* min rectagle query region RQR of a query in query_cell for pivot p */
v_type min_RQR(struct cell * query_cell, unsigned int num_pivots, int p,  v_type dist_threshold)
{
    v_type min_rqr = FLT_MAX;
    long unsigned int num_vectors = 0;

    struct vector * query_vectors = get_sub_cells_vectors_ps(query_cell, num_pivots, &num_vectors);

    for(int i = 0; i < num_vectors; i++)
    {
        // lead vectors stored in cell
        if(dist_threshold - query_vectors[i].values[p] < min_rqr)
                min_rqr = dist_threshold - query_vectors[i].values[p];
    }
    
    for(int i = 0; i < num_vectors; i++)
        free(query_vectors[i].values);
    free(query_vectors);

    return min_rqr;
}
/* check if list of pairs already has query vector */
int has_query_vector(struct pairs * pairs, vector * query_vector)
{
    if(pairs->num_query_vectors == 0)
        return 0;

    for(int i = 0; i < pairs->num_query_vectors; i++)
    {
        // printf("%p vs %p\n", pairs->query_vectors[i], query_vector);
        if
        (
            pairs->query_vectors[i]->table_id == query_vector->table_id
            &&
            pairs->query_vectors[i]->set_id == query_vector->set_id
            &&
            pairs->query_vectors[i]->pos == query_vector->pos
        )
            return i;
        
            
    }

    return -1;
}
/* add candidate pair */
enum response add_candidate_pair(struct pairs * pairs, struct vector * q, struct cell * cell)
{
    if(pairs == NULL)
        exit_with_failure("Error in query_engine.c: NULL pointer to result pairs!");

    unsigned int num_query_vectors = pairs->num_query_vectors;
    int q_idx = has_query_vector(pairs, q);

     // query vector doesn't exist / the first to be added
    if(q_idx == -1 || num_query_vectors == 0)
    {
        if(q_idx == -1)
            q_idx = num_query_vectors; // add at last pos

        // alloc for new query vector
        pairs->query_vectors = realloc(pairs->query_vectors, sizeof(struct vector *) * (num_query_vectors + 1));
        pairs->candidate_pairs = realloc(pairs->candidate_pairs, sizeof(struct candidate_pair *) * (num_query_vectors + 1));
        pairs->matching_pairs = realloc(pairs->matching_pairs, sizeof(struct matching_pair *) * (num_query_vectors + 1));
        
        if(pairs->query_vectors == NULL || pairs->candidate_pairs == NULL) 
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new candidate pair.");
        pairs->query_vectors[q_idx] = q;
        
        pairs->candidate_pairs[q_idx] = malloc(sizeof(struct candidate_pair));
         if(pairs->candidate_pairs[q_idx] == NULL) 
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new candidate pair.");
        
        pairs->candidate_pairs[q_idx]->num_candidates = 0;

        pairs->num_query_vectors++;
    }
    else
    {
        // vector already in pairs, free new duplicate vector
        free(q->values);
        free(q);
    }
    
    // add candidate (todo: check if candidate cell already exists)
    struct candidate_pair * cp = pairs->candidate_pairs[q_idx];
    cp->cells = realloc(cp->cells, sizeof(struct cell *) * (cp->num_candidates + 1));
    if(cp->cells == NULL)
        exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new candidate pair cells.");
    cp->cells[cp->num_candidates] = cell;
    cp->num_candidates++;

    return OK;
}

/* add candidate pair */
enum response add_matching_pair(struct pairs * pairs, struct vector * q, struct cell * cell)
{
    if(pairs == NULL)
        exit_with_failure("Error in query_engine.c: NULL pointer to result pairs!");

    unsigned int num_query_vectors = pairs->num_query_vectors;
    int q_idx = has_query_vector(pairs, q);

    // query vector doesn't exist / the first to be added
    if(q_idx == -1 || num_query_vectors == 0)
    {
        if(q_idx == -1)
            q_idx = num_query_vectors; // add at last pos

        // alloc for new query vector
        pairs->query_vectors = realloc(pairs->query_vectors, sizeof(struct vector *) * (num_query_vectors + 1));
        pairs->candidate_pairs = realloc(pairs->candidate_pairs, sizeof(struct candidate_pair *) * (num_query_vectors + 1));
        pairs->matching_pairs = realloc(pairs->matching_pairs, sizeof(struct matching_pair *) * (num_query_vectors + 1));
        
        if(pairs->query_vectors == NULL || pairs->matching_pairs == NULL) 
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new match pair.");
        pairs->query_vectors[q_idx] = q;

        printf("match pairs %d malloc\n.\n.\n.\n.\n.\n", q_idx);
        pairs->matching_pairs[q_idx] = malloc(sizeof(struct matching_pair));
         if(pairs->matching_pairs[q_idx] == NULL) 
            exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new match pair.");
        
        pairs->matching_pairs[q_idx]->num_match = 0;
        pairs->num_query_vectors++;
    }
    else
    {
        // vector already in pairs, free new duplicate vector
        free(q->values);
        free(q);
    }
    
    
    // add candidate
    struct matching_pair * mp = pairs->matching_pairs[q_idx];
    mp->cells = realloc(mp->cells, sizeof(struct cell *) * (mp->num_match + 1));
    if(mp->cells == NULL)
        exit_with_failure("Error in query_engine.c: Coudn't reallocate memory for new candidate pair cells.");
    mp->cells[mp->num_match] = cell;
    mp->num_match++;

    return OK;
}
/* destroy result pairs */
enum response destroy_pairs(struct pairs * pairs)
{
    if(pairs == NULL)
        exit_with_failure("Error in query_engine.c: NULL pointer to pairs!");

    if(pairs->num_query_vectors == 0)
    {
        free(pairs);
        return OK;
    }

    for(int i = pairs->num_query_vectors - 1; i >= 0; i--)
    {
        if(pairs->candidate_pairs != NULL)
        {
            struct candidate_pair * curr_cp = pairs->candidate_pairs[i];
            free(curr_cp->cells);
            free(pairs->candidate_pairs[i]);
        }
        
        if(pairs->matching_pairs != NULL)
        {
            struct matching_pair * curr_mp = pairs->matching_pairs[i];
            // free(curr_mp->cells);
            free(curr_mp);
        }
        free(pairs->query_vectors[i]->values);
        free(pairs->query_vectors[i]);
    }

    if(pairs->matching_pairs != NULL)
        free(pairs->matching_pairs);
    if(pairs->candidate_pairs != NULL)
        free(pairs->candidate_pairs);
    if(pairs->query_vectors != NULL)
        free(pairs->query_vectors);
        
    free(pairs);
}

struct pairs * init_pairs()
{
    struct pairs * pairs = malloc(sizeof(struct pairs));
    if(pairs == NULL)
        exit_with_failure("Error in pexeso.c: Couldn't allocate memory for pairs!");

    pairs->num_query_vectors = 0u;
    pairs->candidate_pairs = NULL;
    pairs->matching_pairs = NULL;
    pairs->query_vectors = NULL;

    return pairs;
}