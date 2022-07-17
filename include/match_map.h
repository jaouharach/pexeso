#include <stdio.h>

// match and mismatch map
struct match_map {
    struct sid query_set;
    struct sid * sets; // pointer to set in inverted index
    bool * joinable;
    unsigned int * match_count;
    unsigned int * mismatch_count;
    unsigned int * u; // num vectors in column Q that have no matches column in S
    unsigned int total_checked_vectors;
    unsigned int num_sets;
    float query_time;
};

/* create a match/mismatch map for all query sets */
struct match_map * init_match_maps(struct inv_index * index, struct sid * query_sets, int num_query_sets);

/* update match count for a given set */
enum response update_match_count(struct match_map * map_list, int map_idx, struct sid * query_set, int set_idx, float join_threshold, unsigned int query_set_size);

/* update mismatch count for a given set */
enum response update_mismatch_count(struct match_map * map_list, int map_idx, int set_idx);

/* check if set id is in map */
int has_set(struct match_map * map, struct sid * set_id);

/* get idx of match map for a specific query set */
int get_match_map_idx(struct match_map *map, int num_query_sets, struct sid * query_set);

/* print map */
void dump_match_map_to_console(struct match_map * map, unsigned int map_idx);

/* destroy match map */
enum response match_maps_destroy(struct match_map *map, int num_query_sets);
