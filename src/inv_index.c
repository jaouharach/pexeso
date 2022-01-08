#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/cell.h"
#include "../include/inv_index.h"

/* add entry to inverted index */
enum response inv_index_append_entry(struct inv_index * index, unsigned int table_id, unsigned int set_pos, struct cell * cell)
{
    if(index == NULL)
        exit_with_failure("Error in in inv_index.c: Cannot add entry to NULL index.");

    int entry_idx = has_set(index, table_id, set_pos);

    // empty inverted index
    if(index->num_entries == 0 || entry_idx == -1)
    {
        unsigned int num_entries = index->num_entries;
        index->entries = realloc(index->entries, sizeof(struct entry) * (num_entries + 1));
        if(index->entries == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for new entry");

        struct entry * new_entry = &index->entries[num_entries];
        
        new_entry->id = malloc(sizeof(struct sid));
        new_entry->cells = malloc(sizeof(struct cell*));
        if(new_entry->cells == NULL || new_entry->id == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for new entry members");
        

        new_entry->id->table_id = table_id;
        new_entry->id->set_pos = set_pos;

        new_entry->cells[0] = cell;
        new_entry->num_cells = 1;

        index->num_entries++;

        return OK;
    }
    // set already exists in index
    else
    {   
        // if pair set_id --> {cell}
        if(has_entry(index, entry_idx, cell))
            return OK;

        unsigned int num_cells = index->entries[entry_idx].num_cells;
        index->entries[entry_idx].cells = realloc(index->entries[entry_idx].cells, 
                                                    sizeof(struct cell*) * (num_cells + 1));
        if(index->entries[entry_idx].cells == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for entry's new cells");
        
        index->entries[entry_idx].cells[num_cells] = cell;
        index->entries[entry_idx].num_cells++;

        return OK;
    }
    return FAILED;
}

/* check if index has entry with set_id */
int has_set(struct inv_index * index, unsigned int table_id, unsigned int set_pos)
{
    unsigned int num_entries = index->num_entries;
    if(num_entries == 0)
        return -1;

    // search from last to first (assuming that vectors are inserted in a increasing order of set_id)
    for(int e = num_entries - 1; e >= 0; e--)
    {
        struct entry * curr_entry = &index->entries[e];
        if(curr_entry->id->set_pos == set_pos && curr_entry->id->table_id == table_id)
            return e;
    }
    
    return -1;   
}

/* check if index has entry with set_id */
bool has_entry(struct inv_index * index, unsigned int entry_idx, struct cell * cell)
{
    struct entry * entry = &index->entries[entry_idx];
    for(int c = entry->num_cells - 1; c >= 0; c--)
    {
        // both cell and entry->cells[c] point to the same object
        if
        (
            (entry->cells[c] <= cell && entry->cells[c] >= cell) == true
            &&
            (entry->cells[c] < cell) == false 
            &&
            (entry->cells[c] > cell) == false
        )
            return true;
    }
    return false; 
}

/* print inverted index */
enum response dump_inv_index_to_console(struct inv_index *index)
{
    printf("\t\t::  DISPLAY INVERTED INDEX  ::\n");
    printf("\t\t|       |       |       |\n");
    printf("\t\t|       |       |       |\n");
    printf("\t\t|       |       |       |\n");
    printf("\t\tV       V       V       V\n\n\n");
    for(int e = 0; e < index->num_entries; e++)
    {
        struct entry * curr_entry = &index->entries[e];

        printf("%d: (%u, %u) => {", e, curr_entry->id->table_id, curr_entry->id->set_pos);
        for(int c = 0; c < curr_entry->num_cells; c++)
        {
            if(c == curr_entry->num_cells - 1)
            {
                printf("%p", curr_entry->cells[c]);
                break;
            }
            printf("%p :: ", curr_entry->cells[c]);
        }
        printf("}\n");
    }
    printf("\t\t::  END OF HGRID  ::\n");

    return OK;
}

/* destroy inverted index */
enum response inv_index_destroy(struct inv_index *index)
{   
    // destroy entries
    for (int e = index->num_entries - 1; e >= 0; e--)
    {
        // free entry cells
        free(index->entries[e].cells);
        free(index->entries[e].id);
    }
    free(index->entries);
    free(index);

    return OK;
}
