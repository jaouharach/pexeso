#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/map.h"
#include "../include/inv_index.h"

/* create a match and mismatch map for all sets in the inverted index */
enum response init_match_map(struct inv_index * index, struct match_map * map)
{
    map->num_sets = index->num_entries;
    map->id = malloc(sizeof(struct sid *) * index->num_entries); // points to set in inverted index
    map->match_count = calloc(index->num_entries, sizeof(unsigned int));
    map->mismatch_count = calloc(index->num_entries, sizeof(unsigned int));
    
    for(int i = 0; i < index->num_entries; i++)
    {
        // link entry in map with entry in inverted index
        map->id[i] = index->entries[i].id;
    }
    return OK;
}

/* update match count for a given set */
enum response update_match_count(struct match_map * map, struct sid * set_id)
{
    int set_idx = has_sid(map, set_id);
    if(set_idx == -1)
        exit_with_failure("Error in map.c: Couldn't update match map, set id does not exist!");
    map->match_count[set_idx]++;
}

/* update mismatch count for a given set */
enum response update_mismatch_count(struct match_map * map, struct sid * set_id)
{
    int set_idx = has_sid(map, set_id);
    if(set_idx == -1)
        exit_with_failure("Error in map.c: Couldn't update mismatch map, set id does not exist!");
    map->mismatch_count[set_idx]++;
}

/* check if set id is in map */
int has_sid(struct match_map * map, struct sid * id)
{
    if(map->num_sets == 0)
        return -1;

    for(int s = 0; s < map->num_sets; s++)
    {
        if(map->id[s]->table_id == id->table_id && map->id[s]->set_pos == id->set_pos)
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
    for(int s = 0; s < map->num_sets; s++)
    {
        struct sid * curr_set = map->id[s];

        printf("\t%d: (%u, %u) => {", s, curr_set->table_id, curr_set->set_pos);
        printf("match = %u, mistach = %u}\n", map->match_count[s], map->mismatch_count[s]);
       
    }
    printf("\n\t>>>  END OF MAP  <<<\n\n\n");

}

/* destroy match map */
enum response match_map_destroy(struct match_map *map)
{   
    free(map->id);
    free(map->match_count);
    free(map->mismatch_count);
    free(map);
    return OK;
}
