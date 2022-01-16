#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/inv_index.h"
#include "../include/cell.h"

/* add entry to inverted index */
enum response inv_index_append_entry(struct inv_index * index, struct cell * cell, unsigned int table_id, unsigned int set_pos)
{
    if(index == NULL)
        exit_with_failure("Error in in inv_index.c: Cannot add entry to NULL index.");

    int entry_idx = has_cell(index, cell);
    
    // empty inverted index or first time creating entry for cell
    if(index->num_entries == 0)
    {
        unsigned int num_entries = index->num_entries;
        index->entries = malloc(sizeof(struct entry));
        if(index->entries == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for new entry");

        struct entry * new_entry = &index->entries[num_entries];
        
        new_entry->sets = malloc(sizeof(struct sid *));
        index->sets = malloc(sizeof(struct sid));

        if(new_entry->sets == NULL || index->sets == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for first entry");
        
        // add set id to list of distinct sed ids
        index->sets[0].table_id = table_id;
        index->sets[0].set_pos = set_pos;
        index->num_distinct_sets = 1;

        // link set in cell entry to set in list of distinct sed ids
        new_entry->cell = cell;
        new_entry->sets[0] = &index->sets[0];
        new_entry->num_sets = 1;

        index->num_entries++;

        return OK;
    }
    // inverted index is not empty

    // check if set id has not been previously inserted (for another cell)
    int set_idx = previously_indexed_set(index, table_id, set_pos);
    // if set not in inverted index, add set
    if(set_idx == -1)
    {
        unsigned int num_distinct_sets = index->num_distinct_sets;
        // create set in list of distinct sets 
        index->sets = realloc(index->sets, sizeof(struct sid) * (num_distinct_sets + 1));
        if(index->sets == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for entry's new sets");
        
        index->sets[num_distinct_sets].table_id = table_id;
        index->sets[num_distinct_sets].set_pos = set_pos;
        set_idx = num_distinct_sets; 

        index->num_distinct_sets++;
    }

    if (entry_idx == -1) // a new entry (new leaf cell)
    {
        int num_entries = index->num_entries;
        entry_idx = num_entries;
        // allocate memory fe new entry
        index->entries = realloc(index->entries, sizeof(struct entry) * (index->num_entries + 1));
        if(index->entries == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for new cell entry");
        
        struct entry * new_entry = &index->entries[entry_idx];

        new_entry->cell = cell;
        index->entries[entry_idx].sets = malloc(sizeof(struct sid *));
        if(index->entries[entry_idx].sets == NULL)
            exit_with_failure("Error in inv_index.c: could'nt allocate memory for new entry's sets!");
        
        index->entries[entry_idx].sets[0] = &index->sets[set_idx];
        index->entries[entry_idx].num_sets = 1;

        index->num_entries++;

        return OK;
    }
    else // cell has entry in inverted index, append set to cell entry
    {
        // if pair cell --> {set_id} exists don't add set_id to entry
        if(entry_has_set(index, entry_idx, table_id, set_pos))
            return OK;

        unsigned int num_sets = index->entries[entry_idx].num_sets;
        index->entries[entry_idx].sets = realloc(index->entries[entry_idx].sets, 
                                                    sizeof(struct sid *) * (num_sets + 1));
        if(index->entries[entry_idx].sets == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for entry's new sets");

        index->entries[entry_idx].sets[num_sets] = &index->sets[set_idx];
        index->entries[entry_idx].num_sets++;

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
    for(int s = cell_entry->num_sets - 1; s >= 0; s--)
    {
        if ((cell_entry->sets[s]->table_id == table_id) && (cell_entry->sets[s]->set_pos == set_pos))
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
        struct sid * curr_set = &index->sets[s];
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

        printf("\t%d: %p => {", e, curr_entry->cell);
        for(int s = 0; s < curr_entry->num_sets; s++)
        {
            if(s == curr_entry->num_sets - 1)
            {
                printf("(%u, %u)", curr_entry->sets[s]->table_id, curr_entry->sets[s]->set_pos);
                break;
            }
            printf("(%u, %u) :: ", curr_entry->sets[s]->table_id, curr_entry->sets[s]->set_pos);
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
    free(index->sets);
    free(index);

    return OK;
}
