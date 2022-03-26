#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/inv_index.h"
#include "../include/cell.h"

/* add entry to inverted index (cell -> {(table_id, set_pos)}) */
enum response inv_index_append_entry(struct inv_index * index, struct cell * cell, unsigned int table_id, unsigned int set_pos)
{
    if(index == NULL)
        exit_with_failure("Error in in inv_index.c: Cannot add entry to NULL index.");

    int entry_idx = has_cell(index, cell);
    
    // empty inverted index
    if(index->num_entries == 0)
    {
        if(index->num_distinct_sets != 0)
            index->num_distinct_sets = 0;

        // first cell entry
        index->entries = malloc(sizeof(struct entry));
        if(index->entries == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for new entry");

        struct entry * new_entry = &index->entries[0];
        
        // first set
        new_entry->sets = malloc(sizeof(unsigned long long)); // holds idx of set (tableid, setpos) in list of distinct sets
        index->distinct_sets = malloc(sizeof(struct sid));

        if(new_entry->sets == NULL || index->distinct_sets == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for first entry");
        
        // add set to list of distinct set ids
        index->distinct_sets[0].table_id = table_id;
        index->distinct_sets[0].set_pos = set_pos;
        
        // link set in cell entry to set in list of distinct sed ids
        new_entry->cell = cell;
        new_entry->sets[0] = 0; // first cell -> {first set}
        new_entry->num_sets = 1;

        index->num_distinct_sets = 1;
        index->num_entries = 1;

        return OK;
    }

    // inverted index is not empty
    // check if set id has not been previously inserted (for another cell)
    long long int set_idx = previously_indexed_set(index, table_id, set_pos);
    if(set_idx == -1) // if set not in inverted index, add set
    {
        unsigned int num_distinct_sets = index->num_distinct_sets;
        // create set in list of distinct sets 

        index->distinct_sets = realloc(index->distinct_sets, sizeof(struct sid) * (num_distinct_sets + 1));
        if(index->distinct_sets == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for entry's new sets");
        index->distinct_sets[num_distinct_sets].table_id = table_id;
        index->distinct_sets[num_distinct_sets].set_pos = set_pos;

        
        set_idx = index->num_distinct_sets; 
        index->num_distinct_sets++;
    }

    if (entry_idx == -1) // a new entry (new leaf cell)
    {
        int num_entries = index->num_entries;
        entry_idx = num_entries;
        // allocate memory fe new entry
        struct entry * entries = realloc(index->entries, sizeof(struct entry) * (index->num_entries + 1));
        if(entries == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for new cell entry");
        
        index->entries = entries;
        
        struct entry * new_entry = &index->entries[entry_idx];

        new_entry->cell = cell;
        new_entry->sets = malloc(sizeof(unsigned long long));
        if(new_entry->sets == NULL)
            exit_with_failure("Error in inv_index.c: could'nt allocate memory for new entry's sets!");
        
        // first set in entry
        new_entry->sets[0] = set_idx;
        new_entry->num_sets = 1;
        index->num_entries++;        

        return OK;
    }
    else // cell has entry in inverted index, append set to cell entry
    {
        // if pair cell --> {set_id} exists don't add set_id to entry
        if(entry_has_set(index, entry_idx, table_id, set_pos))
            return OK;

        struct entry * curr_entry = &index->entries[entry_idx];

        unsigned int num_sets = curr_entry->num_sets;
        curr_entry->sets = realloc(curr_entry->sets, 
                                                    sizeof(unsigned long long) * (num_sets + 1));
        if(curr_entry->sets == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for entry's new sets");

        curr_entry->sets[curr_entry->num_sets] = set_idx;
        curr_entry->num_sets++;        

        return OK;  
    }
    
    return FAILED;
}

/* check if index has entry for cell  */
int has_cell(struct inv_index * index, struct cell * cell)
{ 
    unsigned int num_entries = index->num_entries;
    if(num_entries == 0)
        return -1;

    for(int e = num_entries - 1; e >= 0; e--)
    {
        struct entry * curr_entry = &index->entries[e];
        // both cell and curr_entry->cell point to the same object
        if(pointer_cmp(curr_entry->cell, cell))
            return e;
    }
    
    return -1;   
}

/* check if cell entry has set_id */
bool entry_has_set(struct inv_index * index, unsigned int entry_idx, unsigned int table_id, unsigned int set_pos)
{
    struct entry * cell_entry = &index->entries[entry_idx];

    if(cell_entry->num_sets == 0)
        return false;

    unsigned long long set_idx;
    for(int s = cell_entry->num_sets - 1; s >= 0; s--)
    {
        set_idx = cell_entry->sets[s];
        if ((index->distinct_sets[set_idx].table_id == table_id) && (index->distinct_sets[set_idx].set_pos == set_pos))
            return true;
    }
    return false; 
}

/* check if set id has already been inserted into a cell entry */
int previously_indexed_set(struct inv_index * index, unsigned int table_id, unsigned int set_pos)
{
    unsigned int num_distinct_sets = index->num_distinct_sets;
    if(num_distinct_sets == 0)
        return -1;

    for(int s = num_distinct_sets - 1; s >= 0; s--)
    {
        struct sid * curr_set = &index->distinct_sets[s];
        // both cell and curr_entry->cell point to the same object
        if(curr_set->table_id == table_id && curr_set->set_pos == set_pos)
            return s;
    }

    return -1;   
}


/* print inverted index */
void dump_inv_index_to_console(struct inv_index *index)
{
    printf("\n\n\n\t...............................\n");
    printf("\t::      INVERTED INDEX       ::\n");
    printf("\t...............................\n\n\n");
    for(int e = 0; e < index->num_entries; e++)
    {
        struct entry * curr_entry = &index->entries[e];
        unsigned long long set_idx;
        printf("\t%d: %p => {", e, curr_entry->cell);
        for(int s = 0; s < curr_entry->num_sets; s++)
        {
            set_idx = curr_entry->sets[s];
            if(s == curr_entry->num_sets - 1)
            {
                printf("(%u, %u)", index->distinct_sets[set_idx].table_id, index->distinct_sets[set_idx].set_pos);
                break;
            }
            printf("(%u, %u) :: ", index->distinct_sets[set_idx].table_id, index->distinct_sets[set_idx].set_pos);
        }
        printf("}\n");
    }
    printf("\n\t>>>  END OF INVERTED INDEX  <<<\n\n\n");

}

/* destroy inverted index */
enum response inv_index_destroy(struct inv_index *index)
{   
    // destroy entries
    for (int e = index->num_entries - 1; e >= 0; e--)
    {
        // free entry sets
        free(index->entries[e].sets);
    }
    free(index->entries);
    free(index->distinct_sets);
    free(index);

    return OK;
}
