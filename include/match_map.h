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
    unsigned int num_dist_calc;
    unsigned int * has_match_for_curr_qvec; // a flag that is used by query vectors to assess whether they have a match in a set or no, since set vectors can be stored in multiple cells and the algorithm performs quering cellwise.
    unsigned int num_sets;
    float query_time;
};

/* create a match/mismatch map for all query sets */
struct match_map * init_match_maps(struct inv_index * index, struct sid * query_sets, int num_query_sets);

/* reset flag query vector has match in curr vector*/
enum response reset_has_match_flag(struct match_map * match_map);

/* update |U|+= 1,  for every set that doesn't have a match for curr query vector*/
enum response update_zero_match_counter(struct match_map * match_map);

/* update match count for a given set */
enum response update_match_count(struct match_map * map_list, int map_idx, struct sid * query_set, int set_idx, float join_threshold, unsigned int query_set_size);

/* update mismatch count for a given set */
enum response update_mismatch_count(struct match_map * map_list, int map_idx, int set_idx);

/* check if set id is in map */
unsigned long get_set_idx(struct match_map * map, struct sid * set_id);

/* get idx of match map for a specific query set */
int get_match_map_idx(struct match_map *map, int num_query_sets, struct sid * query_set);

/* print map */
void dump_match_map_to_console(struct match_map * map, unsigned int map_idx);

/* check if match map has a sorted list of set-ids */
bool is_sorted(struct match_map * map);

/* destroy match map */
enum response match_maps_destroy(struct match_map *map, int num_query_sets);