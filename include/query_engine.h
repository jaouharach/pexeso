#include <stdio.h>
#include "../include/inv_index.h"

// a query result returns all cells with matching vectors to the query vector.
struct query_result
{
    unsigned int num_matches;
    v_type min_distance;
    struct matching_pair * mpair;
    unsigned int *vid_pos;
    struct vid *vector_ids;
};

struct query_settings
{
    unsigned int join_threshold; // T 
    v_type dist_threshold; // tau
};

struct matching_pair
{
    struct vector ** query_vectors; // query vectors
    struct cell ** cells;
    unsigned int num_match; // number of matching cells.
};

struct candidate_pair
{
    struct  vector ** query_vectors; // query vectors
    struct cell ** cells; // pointers to cells in the index
    unsigned int num_candidates; // number of candidate cells (cells that couldn't be filtered)
};

/* get candidate and matching cells of a query cell */
enum response block(struct cell *query_cell, struct cell * r_cell, 
            struct matching_pair * mpair, struct candidate_pair * cpair, struct grid_settings * settings);

enum response verify(struct grid * grid, struct matching_pair * mpair, struct candidate_pair * cpair,
            struct inv_index * index, struct match_map * match_map, unsigned int query_set_size);
            
/* initialize query settings */
struct query_settings * init_query_settings(v_type dist_threshold, unsigned int join_threshold);

/* 
    Given two vectors q and x, a set P of pivot vectors, a distance function d, 
    and a threshold τ.
    if q matches x, then d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
    note: q and v are in pivot space.
    return OK if vector doesn't satisfy: d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
*/
enum response pivot_filter(struct vector * q, struct vector * x, 
                            unsigned int num_pivots, v_type dist_threshold);
                            
/* 
    Given two vectors q and x, a set P of pivot vectors, a distance function d, 
    and a threshold τ, if there exists a pivot p ∈ P such that d(x, p) +d(q, p) ≤ τ,
    then q matches x.
*/
enum response pivot_match(struct vector * q, struct vector * x, 
                            unsigned int num_pivots, v_type dist_threshold);
/* 
   Given a cell c and a mapped query vector q in the pivot space,
   if c ∩ SQR(q, τ) = ∅, then for any mapped vector x ∈ c, 
   its original vector x does not match q.

   function returns OK if cell c can be filtered (pruned).
*/
enum response vector_cell_filter(struct vector * q, struct cell * cell, 
                            unsigned int num_pivots, v_type dist_threshold);
               
/* 
   Given a target cell c and a query cell cq in the pivot space, 
   if c ∩ SQR(cq.center, (τ+cq.length / 2)) = ∅, 
   then for any mapped vector x ∈ c and any query vector q ∈ cq, 
   their original vectors do not match.

   function returns OK if cell c can be filtered (pruned).
*/
enum response cell_cell_filter(struct cell * cell, struct cell * query_cell, 
                            unsigned int num_pivots, v_type dist_threshold);

/* 
    Given a target cell c and a mapped query vector q in the pivot space,
    if there exists a pivot p ∈ P such that c ∩ RQR(q, p, τ) = c, 
    then for any vector x ∈ c, the original vector x matches the query vector q.

   function returns OK if cell c is a match to q.
*/
enum response vector_cell_match(struct vector * q, struct cell * cell, 
                            unsigned int num_pivots, v_type dist_threshold);

/* 
   Given a target cell c and a query cell cq in the pivot space, 
   if there exists a pivot p ∈ P such that c ∩ min(RQR(q, p, τ)) = c, 
   then for any mapped vector x ∈ c and any query vector q ∈ cq,
   their original vectors match.

   function returns OK if cell c is a match to cq.
*/
enum response cell_cell_match(struct cell * cell, struct cell * query_cell, 
                            unsigned int num_pivots, v_type dist_threshold);

/* min rectagle query region RQR of a query in query_cell for pivot p */
v_type min_RQR(struct cell * query_cell, unsigned int num_pivots, int p,  v_type dist_threshold);

/* add candidate pair */
enum response add_candidate_pair(struct candidate_pair *cpair, struct vector * query_vector, struct cell *candidate);

/* add candidate pair */
enum response add_matching_pair(struct matching_pair *mpair, struct vector * query_vector, struct cell *match);

/* destroy candidate pair */
enum response destroy_candidate_pairs(struct candidate_pair *cpair);

/* destroy matching pair */
enum response destroy_matching_pairs(struct matching_pair *mpair);