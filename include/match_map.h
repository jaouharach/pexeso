#include <stdio.h>

// match and mismatch map
struct match_map {
    struct entry ** set_entry; // pointer to set entry in inverted index
    bool * joinable;
    unsigned int * match_count;
    unsigned int * mismatch_count;
    unsigned int num_sets;
};

/* create a match/mismatch map for all sets in the inverted index */
enum response init_match_map(struct inv_index * index, struct match_map * map);

/* update match count for a given set */
enum response update_match_count(struct match_map * map, struct entry * set_entry);

/* update mismatch count for a given set */
enum response update_mismatch_count(struct match_map * map, struct entry * set_entry);

/* check if set id is in map */
int has_set_entry(struct match_map * map, struct entry * set_entry);
/* print map */
void dump_match_map_to_console(struct match_map * map);

/* destroy match map */
enum response match_map_destroy(struct match_map *map);