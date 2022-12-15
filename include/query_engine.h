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
    float join_threshold; // T 
    float dist_threshold; // tau
    int num_query_sets;
    int min_query_set_size;
    int max_query_set_size;
    unsigned int mtr_buffered_memory_size; // amount of memory (in MB) to store metric space vectors for query grid
    unsigned int ps_buffered_memory_size; // amount of memory (in MB) to store pivot space vectors for query grid
};

// pairs of candidate and matching cells for every query vector
/* pairs[i] = <q', has_candidates, has_matches, {candidate_cells}, {matching_cells}>  */
struct pairs {
    struct vector * query_vectors;
    struct vector * query_vectors_mtr;
    struct matching_pair * matching_pairs;
    struct candidate_pair * candidate_pairs;
    bool * has_candidates;
    bool * has_matches;
    unsigned int num_pairs;
};
struct matching_pair
{
    struct cell **cells;
    unsigned int num_match; // number of matching cells.
};

struct candidate_pair
{
    struct cell ** cells; // pointers to cells in the index
    unsigned int num_candidates; // number of candidate cells (cells that couldn't be filtered)
};


/* get candidate and matching cells of a query cell */
enum response block(struct cell *query_cell, struct cell * r_cell, 
                    struct pairs * pairs, struct grid_settings * settings, struct match_map * match_map);

enum response verify(struct grid * grid, struct pairs * pairs,
            struct inv_index * index, struct match_map * match_map);
            
/* initialize query settings */
struct query_settings * init_query_settings(v_type dist_threshold, float join_threshold, int num_query_sets, int min_query_set_size, int max_query_set_size, 
                                            double mtr_buffered_memory_size);

/* check if candidate pair already exists */
int is_candidate_cell(struct pairs * pairs, struct vector * q, struct cell * cell, unsigned int num_pivots);

/* check if a vector is in the SQR of a query vector */
bool vecotr_in_SQR(struct vector * v, struct vector * q, unsigned int num_pivots, float dist_threshold);

/* check if a vector is in the RGR of a query vector and pivot p */
bool vector_in_RQR(struct vector * v, struct vector * q, unsigned int p, float dist_threshold);

/* 
    Given two vectors q and x, a set P of pivot vectors, a distance function d, 
    and a threshold τ.
    if q matches x, then d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
    note: q and v are in pivot space.
    return OK if vector doesn't satisfy: d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ
*/
enum response pivot_filter(struct vector * q, struct vector * x, 
                            unsigned int num_pivots, float dist_threshold);
                            
/* 
    Given two vectors q and x, a set P of pivot vectors, a distance function d, 
    and a threshold τ, if there exists a pivot p ∈ P such that d(x, p) +d(q, p) ≤ τ,
    then q matches x.
*/
enum response pivot_match(struct vector * q, struct vector * x, 
                            unsigned int num_pivots, float dist_threshold);
/* 
   Given a cell c and a mapped query vector q in the pivot space,
   if c ∩ SQR(q, τ) = ∅, then for any mapped vector x ∈ c, 
   its original vector x does not match q.

   function returns OK if cell c can be filtered (pruned).
*/
enum response vector_cell_filter(struct vector * q, struct cell * cell, struct vector * cell_vectors,
                            struct grid_settings * settings, float dist_threshold);
               
/* 
   Given a target cell c and a query cell cq in the pivot space, 
   if c ∩ SQR(cq.center, (τ+cq.length / 2)) = ∅, 
   then for any mapped vector x ∈ c and any query vector q ∈ cq, 
   their original vectors do not match.

   function returns OK if cell c can be filtered (pruned).
*/
enum response cell_cell_filter(struct cell * cell, struct cell * query_cell, 
                            struct grid_settings * settings, float dist_threshold);

/* 
    Given a target cell c and a mapped query vector q in the pivot space,
    if there exists a pivot p ∈ P such that c ∩ RQR(q, p, τ) = c, 
    then for any vector x ∈ c, the original vector x matches the query vector q.

   function returns OK if cell c is a match to q.
*/
enum response vector_cell_match(struct vector * q, struct cell * cell, struct vector * cell_vectors,
                            struct grid_settings * settings, float dist_threshold);

/* 
   Given a target cell c and a query cell cq in the pivot space, 
   if there exists a pivot p ∈ P such that c ∩ min(RQR(q, p, τ)) = c, 
   then for any mapped vector x ∈ c and any query vector q ∈ cq,
   their original vectors match.

   function returns OK if cell c is a match to cq.
*/
enum response cell_cell_match(struct cell * cell, struct cell * query_cell, 
                            struct grid_settings * settings, float dist_threshold);

/* min rectagle query region RQR of a query vector  */
int min_RQR(struct vector *q, unsigned int num_pivots, float dist_threshold);

/* add candidate pair */
enum response add_candidate_pair(struct pairs * pairs, struct vector * query_vector, struct vector * query_vector_mtr, struct cell *candidate, unsigned int num_pivots, unsigned int mtr_vector_length);

/* add candidate pair */
enum response add_matching_pair(struct pairs * pairs, struct vector * query_vector, struct vector * query_vector_mtr, struct cell *match, unsigned int num_pivots, unsigned int mtr_vector_length);

/* check if list of pairs already has query vector */
int has_query_vector(struct pairs * pairs, vector * query_vector);

/* destroy result pairs */
enum response destroy_pairs(struct pairs * pairs);

/* init pairs */
struct pairs * init_pairs();