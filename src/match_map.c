#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/match_map.h"
#include "../include/inv_index.h"

/* reset flag query vector has match in curr vector*/
enum response reset_has_match_flag(struct match_map * match_map)
{
    for(unsigned long s = 0; s < match_map->num_sets; s++)
    {
        match_map->has_match_for_curr_qvec[s] = 0;
    }
    return OK;
}
/* update |U|+= 1,  for every set that doesn't have a match for curr query vector*/
enum response update_zero_match_counter(struct match_map * match_map)
{
    for(unsigned long s = 0; s < match_map->num_sets; s++)
    {
        if(match_map->has_match_for_curr_qvec[s] == 0)
        {
            match_map->u[s] += 1;
        }
    }
    // reset flag for future query vector
    if(reset_has_match_flag(match_map))
        return OK;
    else 
        return FAILED;
}
/* create a match/mismatch map for all query sets */
struct match_map * init_match_maps(struct inv_index * index, struct sid * query_sets, int num_query_sets)
{
    if(num_query_sets < 0)
        exit_with_failure("Error in match_map.c: cannot initialize match map for zero query sets!");

    struct match_map * match_map = malloc(sizeof(struct match_map) * num_query_sets);
    if(match_map == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for match maps.");

    for(int i = 0; i < num_query_sets; i++)
    {
        struct match_map * curr_map = &match_map[i];
        curr_map->query_set.table_id = query_sets[i].table_id;
        curr_map->query_set.set_id = query_sets[i].set_id;
        curr_map->query_set.set_size = query_sets[i].set_size;

        curr_map->num_sets = index->num_distinct_sets;
        curr_map->sets = (struct sid *)malloc(sizeof(struct sid) * index->num_distinct_sets); // points to set in inverted index
        curr_map->match_count = calloc(index->num_distinct_sets, sizeof(unsigned int));
        curr_map->mismatch_count = calloc(index->num_distinct_sets, sizeof(unsigned int));
        curr_map->u = calloc(index->num_distinct_sets, sizeof(unsigned int));
        curr_map->has_match_for_curr_qvec = calloc(index->num_distinct_sets, sizeof(unsigned int));
        curr_map->joinable = calloc(index->num_distinct_sets, sizeof(bool));
        curr_map->total_checked_vectors = 0;
        curr_map->num_dist_calc = 0;
        curr_map->query_time = 0;

        for(int s = 0; s < index->num_distinct_sets; s++)
        {
            struct sid * curr_set = &curr_map->sets[s];
            // from inverted index copy set ids to map  
            if(index->distinct_sets == NULL)
                exit_with_failure("Error in match_map.c: NULL pointer! couldn't link entry in matchmap to entry in inverted index!");
            
            curr_set->table_id = index->distinct_sets[s].table_id;
            curr_set->set_id = index->distinct_sets[s].set_id;
            curr_set->set_size = index->distinct_sets[s].set_size;
            curr_map->match_count[s] = 0;
            curr_map->mismatch_count[s] = 0;
            curr_map->u[s] = 0;
            curr_map->has_match_for_curr_qvec[s] = 0;
            curr_map->joinable[s] = false;
        }
    }
    
    // free memory
    free(query_sets);
    return match_map;
}

/* update match count for a given set */
enum response update_match_count(struct match_map * map_list, int map_idx, struct sid * query_set, int set_idx, float join_threshold, unsigned int query_set_size)
{
    map_list[map_idx].match_count[set_idx] = map_list[map_idx].match_count[set_idx] + 1;

    if(map_list[map_idx].match_count[set_idx] >= ceil(join_threshold * query_set_size))
    {
        map_list[map_idx].joinable[set_idx] = true;
    }
    return OK;
}

/* update mismatch count for a given set */
enum response update_mismatch_count(struct match_map * map_list, int map_idx, int set_idx)
{
    map_list[map_idx].mismatch_count[set_idx] = map_list[map_idx].mismatch_count[set_idx] + 1;
    return OK;
}

/* check if set id is in map */
// unsigned long get_set_idx(struct match_map * map, struct sid * sid)
// {
//     int bs_result = -2;
//     int seq_result = -3;

//     if(map->num_sets == 0)
//         return -1;

//     // first occurence of sid->table_id
//     bs_result =  binary_search(map->sets, sid, 0, map->num_sets - 1);

//     // perform sequential search
//     for(int s = 0; s < map->num_sets; s++)
//     {
//         if(map->sets[s].table_id == sid->table_id && map->sets[s].set_id == sid->set_id)
//         {
//             seq_result =  s;
//         }
//     }
    
//     if(bs_result != seq_result)
//     {
//         printf("\nbs result = %d, seq result = %d\n", bs_result, seq_result);
//         if(bs_result < 0)
//             printf("bs found nothing\n");
//         else
//             printf("bs set is (%u, %u)\n", map->sets[bs_result].table_id, map->sets[bs_result].set_id);
//         if(seq_result < 0)
//             printf("seq found nothing\n");
//         else
//             printf("seq set is (%u, %u)\n", map->sets[seq_result].table_id, map->sets[seq_result].set_id);  
//     } 
// }

/* get idx of match map for a specific query set */
int get_match_map_idx(struct match_map *map, int num_query_sets, struct sid * sid)
{
    for(int m = 0; m < num_query_sets; m++)
    {
        if(map[m].query_set.table_id == sid->table_id && map[m].query_set.set_id == sid->set_id)
            return m;
    }
    return -1;
}

/* print map */
void dump_match_map_to_console(struct match_map * map, unsigned int map_idx)
{
    printf("\n\n\n\t............................................................\n");
    printf("\t::  MATCH  & MISMATCH MAP FOR QUERY SET Q: (%u, %u), |Q| = %u ::\n", map[map_idx].query_set.table_id, map[map_idx].query_set.set_id, map[map_idx].query_set.set_size);
    printf("\t \t\t(Query time = %.2f seconds)\t\t\n", map[map_idx].query_time / 1000000);

    printf("\t............................................................\n\n\n");
    struct sid curr_set;
    for(unsigned long s = 0; s < map[map_idx].num_sets; s++)
    {
        if(map[map_idx].joinable[s]) // only print joinable sets
        {
            curr_set = map[map_idx].sets[s];
            printf("\t%ld: \033[1;33m S: (%u, %u) |S| = %u\033[0m => {", s, curr_set.table_id, curr_set.set_id, curr_set.set_size);
            printf("match = %u, mismatch = %u, |U| = %u, joinable: %s}\n", map[map_idx].match_count[s], map[map_idx].mismatch_count[s], map[map_idx].u[s], map[map_idx].joinable[s] ? "\033[0;32mtrue\033[0m" : "\033[1;31mfalse\033[0m");
        }
    }
    printf("\n\t>>>  END OF MAP  <<<\n\n\n");

}

void dump_csv_results_to_console(struct match_map * map, unsigned int map_idx)
{
    printf("\t \t\tQ = (%u, %u), -Q- = %u (Query time = %.2f seconds)\t\t\n", 
    map[map_idx].query_set.table_id, map[map_idx].query_set.set_id, map[map_idx].query_set.set_size, map[map_idx].query_time / 1000000);
    printf("\t............................................................\n\n\n");
    struct sid curr_set;
    printf("tqq,tss,nooverlap\n");
    for(unsigned long s = 0; s < map[map_idx].num_sets; s++)
    {
        if(map[map_idx].joinable[s]) // only print joinable sets
        {
            curr_set = map[map_idx].sets[s];
            printf("t%uc%u,t%uc%u,%u\n", map[map_idx].query_set.table_id, map[map_idx].query_set.set_id, curr_set.table_id, curr_set.set_id, map[map_idx].u[s]);
        }
    }
    printf("\t............................................................\n\n\n");
}


/* destroy match map */
enum response match_maps_destroy(struct match_map *map, int num_query_sets)
{   
    for(int i = 0; i < num_query_sets; i++)
    {
        free(map[i].sets);
        free(map[i].joinable);
        free(map[i].match_count);
        free(map[i].mismatch_count);
        free(map[i].u);
        free(map[i].has_match_for_curr_qvec);
    }
    free(map);
    return OK;
}

