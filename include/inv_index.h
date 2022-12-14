#include <stdio.h>

// a list of set_id -->{c1, c5, 34, ...} sorted by set_id
struct inv_index {
    unsigned int num_entries; // number of leaf cells with an entry in inverted index
    unsigned long num_distinct_sets; // number of unique set ids in inverted index
    struct entry * entries; // cell entries cell --> {set1, set2, ...}
    struct sid * distinct_sets; // unique sets indexed in inverted index
    
};

struct entry {
    struct cell * cell; // one cell
    unsigned long  * sets; // multiple sets, eatch set is the idx of set (column) in inverted index
    unsigned long * vector_count; // number of vectors fro set s in cell c
    unsigned long num_sets;
};

/* add entry to inverted index */
long inv_index_append_entry(struct inv_index * index, struct cell * cell, unsigned int table_id, unsigned int set_id, unsigned int set_size, unsigned int num_pivots);

/* check if index has entry for cell  */
int has_cell(struct inv_index * index, struct cell * cell, unsigned int num_pivots);
/* check if cell entry has set_id */

unsigned long entry_has_set(struct inv_index * index, int entry_idx, unsigned long set_idx);

/* check if set id has already been inserted into a cell entry */
long previously_indexed_set(struct inv_index * index, unsigned int table_id, unsigned int set_id);

/* print inverted index */
void dump_inv_index_to_console(struct inv_index *index);

/* destroy inverted index */
enum response inv_index_destroy(struct inv_index *index);

/* write inv index to disk */
enum response index_write(struct inv_index * index, const char * work_dir);

/* read inv index from disk */
struct inv_index * index_read(const char * work_dir, 
                    struct grid_settings * settings,
                    struct level * root);