#include <stdio.h>

// match and mismatch map
struct match_map {
    struct sid * sets; // pointer to set in inverted index
    bool * joinable;
    unsigned int * match_count;
    unsigned int * mismatch_count;
    unsigned int num_sets;
};

/* create a match/mismatch map for all sets in the inverted index */
enum response init_match_map(struct inv_index * index, struct match_map * map);

/* update match count for a given set */
enum response update_match_count(struct match_map * map, struct sid * set_id);

/* update mismatch count for a given set */
enum response update_mismatch_count(struct match_map * map, struct sid * set_id);

/* check if set id is in map */
int has_set(struct match_map * map, struct sid * set_id);

/* print map */
void dump_match_map_to_console(struct match_map * map);

/* destroy match map */
enum response match_map_destroy(struct match_map *map);