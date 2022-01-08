#include <stdio.h>

// a list of set_id -->{c1, c5, 34, ...} sorted by set_id
struct inv_index {
    unsigned long long num_entries;
    struct entry * entries;
};

struct entry {
    struct sid * id;
    struct cell ** cells;
    unsigned int num_cells;
    // struct vector * vectors;
};

/* add entry to inverted index */
enum response inv_index_append_entry(struct inv_index * index, unsigned int table_id, unsigned int set_pos, struct cell * cell);

/* check if index has entry with set_id and cell*/
bool has_entry(struct inv_index * index, unsigned int entry_idx, struct cell * cell);

/* check if index has entry with set_id */
int has_set(struct inv_index * index, unsigned int table_id, unsigned int set_pos);

/* print inverted index */
enum response dump_inv_index_to_console(struct inv_index *index);

/* destroy inverted index */
enum response inv_index_destroy(struct inv_index *index);