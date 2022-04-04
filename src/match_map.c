#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/match_map.h"
#include "../include/inv_index.h"

/* create a match and mismatch map for all sets in the inverted index */
enum response init_match_map(struct inv_index * index, struct match_map * map)
{
    map->num_sets = index->num_distinct_sets;
    map->sets = (struct sid *)malloc(sizeof(struct sid) * index->num_distinct_sets); // points to set in inverted index
    map->match_count = calloc(index->num_distinct_sets, sizeof(unsigned int));
    map->mismatch_count = calloc(index->num_distinct_sets, sizeof(unsigned int));
    map->joinable = calloc(index->num_distinct_sets, sizeof(bool));

    for(int s = 0; s < index->num_distinct_sets; s++)
    {
        // from inverted index copy set ids to map  
        if(index->distinct_sets == NULL)
            exit_with_failure("Error in match_map.c: NULL pointer! couldn't link entry in matchmap to entry in inverted index!");
        
        map->sets[s].table_id = index->distinct_sets[s].table_id;
        map->sets[s].set_id = index->distinct_sets[s].set_id;
        map->sets[s].set_size = index->distinct_sets[s].set_size;
        map->match_count[s] = 0;
        map->mismatch_count[s] = 0;
        map->joinable[s] = false;

    }
    return OK;
}

/* update match count for a given set */
enum response update_match_count(struct match_map * map, struct sid * sid, float join_threshold, unsigned int query_set_size)
{
    int set_idx = has_set(map, sid);
    if(set_idx == -1)
        exit_with_failure("Error in match_map.c: Couldn't update match map, set id does not exist in match map!");
    
    map->match_count[set_idx] = map->match_count[set_idx] + 1;

    if(map->match_count[set_idx] >= ceil(join_threshold * query_set_size))
    {
        // printf("\033[1;31mjoin threshold = %f * %u = %f, match count = %u\033[0m", join_threshold, query_set_size, ceil(join_threshold * query_set_size), map->match_count[set_idx]);
        map->joinable[set_idx] = true;
    }
    return OK;
}

/* update mismatch count for a given set */
enum response update_mismatch_count(struct match_map * map, struct sid * sid)
{
    int set_idx = has_set(map, sid);
    if(set_idx == -1)
        exit_with_failure("Error in match_map.c: Couldn't update mismatch map, set id does not exist in match map!");
    map->mismatch_count[set_idx] = map->mismatch_count[set_idx] + 1;
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

/* print map */
void dump_match_map_to_console(struct match_map * map)
{
    printf("\n\n\n\t..............................\n");
    printf("\t::  MATCH  & MISMATCH MAP   ::\n");
    printf("\t..............................\n\n\n");
    struct sid curr_set;
    for(int s = 0; s < map->num_sets; s++)
    {
        curr_set = map->sets[s];
        printf("\t%d: \033[1;33m S: (%u, %u) |S| = %u\033[0m => {", s, curr_set.table_id, curr_set.set_id, curr_set.set_size);
        printf("match = %u, mismatch = %u, joinable: %s}\n", map->match_count[s], map->mismatch_count[s], map->joinable[s] ? "\033[0;32mtrue\033[0m" : "\033[1;31mfalse\033[0m");
       
    }
    printf("\n\t>>>  END OF MAP  <<<\n\n\n");

}

/* destroy match map */
enum response match_map_destroy(struct match_map *map)
{   
    free(map->sets);
    free(map->joinable);
    free(map->match_count);
    free(map->mismatch_count);
    free(map);
    return OK;
}
