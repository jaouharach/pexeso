#include <stdio.h>
#include <float.h>
#include "../include/index.h"
#include "../include/cell.h"
#include "../include/query_engine.h"

/* 
    Given two vectors q and x, a set P of pivot vectors, a distance function d, 
    and a threshold τ, if there exists a pivot p ∈ P such that d(x, p) +d(q, p) ≤ τ,
    then q matches x.
*/
enum response pivot_match(struct vector * q, struct vector * x, 
                            unsigned int num_pivots, v_type dist_threshold)
{
    // find pivot such that d(q, p) − τ ≤ d(x, p) ≤ d(q, p) + τ

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
        min_rqr = min_RQR(query_cell->center, num_pivots, p, dist_threshold);
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
    struct vector * query_vectors = get_vectors(query_cell, num_pivots);
    for(int i = 0; i < query_cell->cell_size; i++)
    {
        // lead vectors stored in cell
        if(dist_threshold - query_vectors[i].values[p] < min_rqr)
                min_rqr = dist_threshold - query_vectors[i].values[p];
    }
    
    for(int i = 0; i < query_cell->cell_size; i++)
        free(query_vectors[i].values);
    free(query_vectors);

    return min_rqr;
}