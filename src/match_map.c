#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/match_map.h"
#include "../include/inv_index.h"

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
        match_map[i].query_set.table_id = query_sets[i].table_id;
        match_map[i].query_set.set_id = query_sets[i].set_id;
        match_map[i].query_set.set_size = query_sets[i].set_size;

        match_map[i].num_sets = index->num_distinct_sets;
        match_map[i].sets = (struct sid *)malloc(sizeof(struct sid) * index->num_distinct_sets); // points to set in inverted index
        match_map[i].match_count = calloc(index->num_distinct_sets, sizeof(unsigned int));
        match_map[i].mismatch_count = calloc(index->num_distinct_sets, sizeof(unsigned int));
        match_map[i].joinable = calloc(index->num_distinct_sets, sizeof(bool));
        match_map[i].total_checked_vectors = 0;
        match_map[i].query_time = 0;

        for(int s = 0; s < index->num_distinct_sets; s++)
        {
            // from inverted index copy set ids to map  
            if(index->distinct_sets == NULL)
                exit_with_failure("Error in match_map.c: NULL pointer! couldn't link entry in matchmap to entry in inverted index!");
            
            match_map[i].sets[s].table_id = index->distinct_sets[s].table_id;
            match_map[i].sets[s].set_id = index->distinct_sets[s].set_id;
            match_map[i].sets[s].set_size = index->distinct_sets[s].set_size;
            match_map[i].match_count[s] = 0;
            match_map[i].mismatch_count[s] = 0;
            match_map[i].joinable[s] = false;
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
int has_set(struct match_map * map, struct sid * sid)
{
    if(map->num_sets == 0)
        return -1;

    for(int s = 0; s < map->num_sets; s++)
    {
        if(map->sets[s].table_id == sid->table_id && map->sets[s].set_id == sid->set_id)
            return s;
    }

    return -1;
}

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
    for(int s = 0; s < map[map_idx].num_sets; s++)
    {
        if(map[map_idx].joinable[s]) // only print joinable sets
        {
            curr_set = map[map_idx].sets[s];
            printf("\t%d: \033[1;33m S: (%u, %u) |S| = %u\033[0m => {", s, curr_set.table_id, curr_set.set_id, curr_set.set_size);
            printf("match = %u, mismatch = %u, joinable: %s}\n", map[map_idx].match_count[s], map[map_idx].mismatch_count[s], map[map_idx].joinable[s] ? "\033[0;32mtrue\033[0m" : "\033[1;31mfalse\033[0m");
        }
    }
    printf("\n\t>>>  END OF MAP  <<<\n\n\n");

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
    }
    free(map);
    return OK;
}
