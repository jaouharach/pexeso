#include <stdio.h>

// a list of set_id -->{c1, c5, 34, ...} sorted by set_id
struct inv_index {
    unsigned long long num_entries; // number of leaf cells with an entry in inverted index
    unsigned long long num_distinct_sets; // number of unique set ids in inverted index
    struct entry * entries; // cell entries cell --> {set1, set2, ...}
    struct sid * distinct_sets; // unique sets indexed in inverted index
};

// struct entry {
//     struct sid * id;
//     struct cell ** cells;
//     unsigned int num_cells;
//     // struct vector * vectors;
// };

// (todo) reverse index order cell -> sets

struct entry {
    struct cell * cell; // one cell
    unsigned long long  * sets; // multiple sets, eatch set is the idx of set id in inverted index
    unsigned int num_sets;
};

/* add entry to inverted index */
enum response inv_index_append_entry(struct inv_index * index, struct cell * cell, unsigned int table_id, unsigned int set_pos);

/* check if index has entry for cell  */
int has_cell(struct inv_index * index, struct cell * cell);
/* check if cell entry has set_id */

bool entry_has_set(struct inv_index * index, unsigned int entry_idx, unsigned int table_id, unsigned int set_pos);

/* check if set id has already been inserted into a cell entry */
int previously_indexed_set(struct inv_index * index, unsigned int table_id, unsigned int set_pos);

/* print inverted index */
void dump_inv_index_to_console(struct inv_index *index);

/* destroy inverted index */
enum response inv_index_destroy(struct inv_index *index);