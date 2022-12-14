#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/inv_index.h"
#include "../include/cell.h"
#include "../include/stats.h"
#include "../include/level.h"

/* add entry to inverted index (cell -> {(table_id, set_id)}) */
long inv_index_append_entry(struct inv_index * index, struct cell * cell, unsigned int table_id, unsigned int set_id, unsigned int set_size, unsigned int num_pivots)
{
    // printf("New vector in set (%u, %u)\n", table_id, set_id);
    if(index == NULL)
        exit_with_failure("Error in in inv_index.c: Cannot add entry to NULL index.");
    
    // int entry_idx = has_cell(index, cell, num_pivots);
    int entry_idx = cell->index_entry_pos;

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
        new_entry->sets = malloc(sizeof(unsigned long)); // holds idx of set (tableid, setpos) in list of distinct sets
        new_entry->vector_count = malloc(sizeof(unsigned long)); // count vectors of set s in cell c
        index->distinct_sets = malloc(sizeof(struct sid));

        if(new_entry->sets == NULL || index->distinct_sets == NULL || new_entry->vector_count == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for first entry");
        
        // add set to list of distinct set ids
        index->distinct_sets[0].table_id = table_id;
        index->distinct_sets[0].set_id = set_id;
        index->distinct_sets[0].set_size = set_size;

        // link set in cell entry to set in list of distinct sed ids
        new_entry->cell = cell;
        cell->index_entry_pos = 0;
        new_entry->sets[0] = 0; // first cell -> {first set}
        new_entry->vector_count[0] = 1;
        new_entry->num_sets = 1;

        index->num_distinct_sets = 1;
        index->num_entries = 1;

        return 0;
    }

    // inverted index is not empty
    // check if set id has not been previously inserted (for another cell)
    long set_in_idx = previously_indexed_set(index, table_id, set_id);
    // printf("Set (%u, %u) in index ? answer = %ld\n", table_id, set_id, set_in_idx);
    unsigned long set_idx = (unsigned long) set_in_idx;
    bool new_set = false;
    if(set_in_idx == -1) // if set not in inverted index, add set
    {
        new_set = true;
        unsigned long num_distinct_sets = index->num_distinct_sets;
        // create set in list of distinct sets 
        index->distinct_sets = realloc(index->distinct_sets, sizeof(struct sid) * (num_distinct_sets + 1));
        if(index->distinct_sets == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for entry's new sets");
        index->distinct_sets[num_distinct_sets].table_id = table_id;
        index->distinct_sets[num_distinct_sets].set_id = set_id;
        index->distinct_sets[num_distinct_sets].set_size = set_size;

        set_idx = index->num_distinct_sets; 
        index->num_distinct_sets++;
    }
    
    if (entry_idx == -1) // a new entry (new leaf cell)
    {
        unsigned int num_entries = index->num_entries;
        entry_idx = (int) num_entries;
        // allocate memory fe new entry
        struct entry * entries = realloc(index->entries, sizeof(struct entry) * (index->num_entries + 1));
        if(entries == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for new cell entry");
        index->entries = entries;
        index->num_entries++;  

        struct entry * new_entry = &index->entries[entry_idx];
        new_entry->cell = cell;
        
        new_entry->sets = malloc(sizeof(unsigned long));
        new_entry->vector_count = malloc(sizeof(unsigned long));
        if(new_entry->sets == NULL || new_entry->vector_count == NULL)
            exit_with_failure("Error in inv_index.c: could'nt allocate memory for new entry's sets!");
        
        // first set in entry
        new_entry->sets[0] = set_idx;
        new_entry->vector_count[0] = 1;
        new_entry->num_sets = 1;
        cell->index_entry_pos = entry_idx;

        return (long) set_idx;
    }
    else // cell has entry in inverted index, append set to cell entry
    {
        // if pair cell --> {set_id} exists don't add set_id to entry and increase vector count 
        if(!new_set)
        {
            unsigned long s = entry_has_set(index, entry_idx, set_idx);
            if(s != -1)
            {
                struct entry * curr_entry = &index->entries[entry_idx];
                curr_entry->vector_count[s] += 1;
                return (long) set_idx;
            }
        }
        
        // if pair cell --> {set_id} exists add set_id to entry
        struct entry * curr_entry = &index->entries[entry_idx];
        unsigned long num_sets = curr_entry->num_sets;
        curr_entry->sets = realloc(curr_entry->sets, sizeof(unsigned long) * (num_sets + 1));
        curr_entry->vector_count = realloc(curr_entry->vector_count, sizeof(unsigned long) * (num_sets + 1));
        if(curr_entry->sets == NULL || curr_entry->vector_count == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for entry's new sets");

        curr_entry->sets[curr_entry->num_sets] = set_idx;
        curr_entry->vector_count[curr_entry->num_sets] = 1;
        curr_entry->num_sets++;        

        return (long) set_idx;  
    }
    
    return -1;
}

/* check if index has entry for cell  */
int has_cell(struct inv_index * index, struct cell * cell, unsigned int num_pivots)
{ 
    unsigned int num_entries = index->num_entries;
    if(num_entries == 0)
        return -1;

    for(unsigned int e = (num_entries - 1); e >= 0; e--)
    {
        struct entry * curr_entry = &index->entries[e];
        int match = 1;
        // both cell and curr_entry->cell have the same celter
        for(int p = 0; p < num_pivots; p++)
        {
            // for two cells to be the same they must have the same center vector
            if(cell->center->values[p] != curr_entry->cell->center->values[p])
            {
                match = 0;
                break;
            }
        }
        if(match == 1)
            return e;
    }
    
    return -1;   
}

/* check if cell entry has set_id */
unsigned long entry_has_set(struct inv_index * index, int entry_idx, unsigned long set_idx)
{
    if(entry_idx == -1)
        exit_with_failure("Error in inv_index.c: Entry index cannot be -1!");

    struct entry * cell_entry = &index->entries[entry_idx];

    if(cell_entry->num_sets == 0)
        return false;

    for(unsigned long s = cell_entry->num_sets - 1; s >= 0; s--)
    {
        if (cell_entry->sets[s] == set_idx)
            return s;    
        if(s == 0) break;
    }
    return -1; 
}

/* check if set id has already been inserted into a cell entry */
long previously_indexed_set(struct inv_index * index, unsigned int table_id, unsigned int set_id)
{
    unsigned long num_distinct_sets = index->num_distinct_sets;
    if(num_distinct_sets == 0)
        return -1;
    // start from end of array for efficiency
    for(unsigned long s = num_distinct_sets - 1; s >= 0; s--)
    {
        struct sid * curr_set = &index->distinct_sets[s];
        // both cell and curr_entry->cell point to the same object
        if(curr_set->table_id == table_id && curr_set->set_id == set_id)
            return  (long) s;
        if(s == 0) break;
    }

    return -1;   
}

/* print inverted index */
void dump_inv_index_to_console(struct inv_index *index)
{
    printf("\n\n\n\t...............................\n");
    printf("\t::      INVERTED INDEX       ::\n");
    printf("\t...............................\n\n\n");
    for(unsigned int e = 0; e < index->num_entries; e++)
    {
        struct entry * curr_entry = &index->entries[e];
        unsigned long set_idx;
        printf("\t%u: cell %d at %p => {", e, curr_entry->cell->id, curr_entry->cell);
        for(unsigned long s = 0; s < curr_entry->num_sets; s++)
        {
            set_idx = curr_entry->sets[s];
            if(s == curr_entry->num_sets - 1)
            {
                printf("(%u, %u)", index->distinct_sets[set_idx].table_id, index->distinct_sets[set_idx].set_id);
                break;
            }
            printf("(%u, %u) :: ", index->distinct_sets[set_idx].table_id, index->distinct_sets[set_idx].set_id);
        }
        printf("}\n");
    }
    printf("\n\t>>>  END OF INVERTED INDEX  <<<\n\n\n");

}

/* destroy inverted index */
enum response inv_index_destroy(struct inv_index *index)
{   
    // destroy entries
    for (int e = (int)(index->num_entries - 1); e >= 0; e--)
    {
        // free entry sets
        free(index->entries[e].sets);
        free(index->entries[e].vector_count);
    }
    free(index->entries);
    free(index->distinct_sets);
    free(index);

    return OK;
}

/* write inv index to disk */
enum response index_write(struct inv_index * index, const char * work_dir)
{
    /* write inverted index to disk */ 
    // make inverted.idx file
    char *inverted_idx_filename = malloc(sizeof(char) * (strlen(work_dir) + 13));
    if (inverted_idx_filename == NULL)
        exit_with_failure("Error in inv_index.c: Couldn't allocate memory for root file name.");
    inverted_idx_filename = strcpy(inverted_idx_filename, work_dir);
    inverted_idx_filename = strcat(inverted_idx_filename, "inverted.idx\0");

    COUNT_PARTIAL_OUTPUT_TIME_START
    FILE *idx_file = fopen(inverted_idx_filename, "wb");
    COUNT_PARTIAL_OUTPUT_TIME_END
    if (idx_file == NULL)
        exit_with_failure("Error in inv_index.c: Couldn't open inverted index file 'inverted.idx'.");

    unsigned int num_entries = index->num_entries;
    unsigned long num_distinct_sets = index->num_distinct_sets;
    struct entry * entries = index->entries;
    struct sid * distinct_sets = index->distinct_sets;

    COUNT_PARTIAL_OUTPUT_TIME_START
    fwrite(&num_entries, sizeof(unsigned int), 1, idx_file);
    fwrite(&num_distinct_sets, sizeof(unsigned long), 1, idx_file);
    fwrite(distinct_sets, sizeof(struct sid), num_distinct_sets, idx_file);

    
    for(unsigned int e = 0; e < num_entries; e++)
    {
        if(entries[e].num_sets > num_distinct_sets)
        {
            printf("Error in inv_index.c: Number of sets in entry %u = %ld cannot exceed #distinct_sets = %ld in inverted index."
            , e, entries[e].num_sets, num_distinct_sets);
            exit(1);
        }
        int filename_size = strlen(entries[e].cell->filename);
        fwrite(&(entries[e].cell->id), sizeof(unsigned int), 1, idx_file);
        fwrite(entries[e].cell->filename, sizeof(char), filename_size, idx_file);
        fwrite(&(entries[e].cell->cell_size), sizeof(unsigned int), 1, idx_file);
        fwrite(&(entries[e].num_sets), sizeof(unsigned long), 1, idx_file);
        fwrite(entries[e].sets, sizeof(unsigned long), entries[e].num_sets, idx_file);
        fwrite(entries[e].vector_count, sizeof(unsigned long), entries[e].num_sets, idx_file);
    }
    COUNT_PARTIAL_OUTPUT_TIME_END

    free(inverted_idx_filename);
    COUNT_PARTIAL_OUTPUT_TIME_START
    fclose(idx_file);
    COUNT_PARTIAL_OUTPUT_TIME_END

    return OK;
}

/* write inv index to disk */
struct inv_index * index_read(const char * work_dir, 
                    struct grid_settings * settings,
                    struct level * root)
{
    struct level * leaf_level = root;
    while(!leaf_level->is_leaf)
        leaf_level = leaf_level->next;
    
    printf("\n\nReading inverted index... ");
    struct inv_index * index = malloc(sizeof(struct inv_index));
    /* write inverted index to disk */ 
    // make inverted.idx file
    char *inverted_idx_filename = malloc(sizeof(char) * (strlen(work_dir) + 13));
    if (inverted_idx_filename == NULL)
        exit_with_failure("Error in inv_index.c: Couldn't allocate memory for root file name.");
    inverted_idx_filename = strcpy(inverted_idx_filename, work_dir);
    inverted_idx_filename = strcat(inverted_idx_filename, "inverted.idx\0");

    COUNT_PARTIAL_INPUT_TIME_START
    FILE *idx_file = fopen(inverted_idx_filename, "rb");
    COUNT_PARTIAL_INPUT_TIME_END
    if (idx_file == NULL)
        exit_with_failure("Error in inv_index.c: Couldn't open inverted index file 'inverted.idx'.");

    COUNT_PARTIAL_INPUT_TIME_START
    fread(&(index->num_entries), sizeof(unsigned int), 1, idx_file);    
    fread(&(index->num_distinct_sets), sizeof(unsigned long), 1, idx_file);
    COUNT_PARTIAL_INPUT_TIME_END
    index->distinct_sets = malloc(sizeof(struct sid) * index->num_distinct_sets);
    index->entries = malloc(sizeof(struct entry) * index->num_entries);

    if (index->entries == NULL  || index->distinct_sets == NULL)
        exit_with_failure("Error in inv_index.c: Couldn't allocate memory for inverted index (entries/sets).");
    
    COUNT_PARTIAL_INPUT_TIME_START
    fread(index->distinct_sets, sizeof(struct sid), index->num_distinct_sets, idx_file);
    printf("num distinct sets = %ld\n", index->num_distinct_sets);
    COUNT_PARTIAL_INPUT_TIME_END
    
    struct entry * entries = index->entries;
    unsigned int cell_id, cell_size;
    char cell_filename[255] = "";

    
    for(unsigned int e = 0; e < index->num_entries; e++)
    {
        COUNT_PARTIAL_INPUT_TIME_START
        fread(&cell_id, sizeof(unsigned int), 1, idx_file);
        COUNT_PARTIAL_INPUT_TIME_END
        for(int c = 0; c < leaf_level->num_cells; c++)
        {
            if(leaf_level->cells[c].id == cell_id)
            {
                entries[e].cell = &(leaf_level->cells[c]); // link cell in inv index with cell in grid
                leaf_level->cells[c].index_entry_pos = e;
            }
        }
            
        int filename_size = strlen(entries[e].cell->filename);
        COUNT_PARTIAL_INPUT_TIME_START
        fread(cell_filename, sizeof(char), filename_size, idx_file);
        cell_filename[filename_size] = '\0';
        fread(&cell_size, sizeof(unsigned int), 1, idx_file);
        
        if(strcmp(cell_filename, entries[e].cell->filename) != 0 || cell_size != entries[e].cell->cell_size)
        {
            printf("Error in inv_index.c: Expected cell %s but read cell %s.", entries[e].cell->filename, cell_filename);
            exit(1); 
        }
        
        fread(&(entries[e].num_sets), sizeof(unsigned long), 1, idx_file);
        if(entries[e].num_sets > index->num_distinct_sets)
            exit_with_failure("Error in inv_index.c: Number of sets in entry cannot exceed #distinct sets in inverted index.");

        // printf("cell %d, entry pos = %d, num sets = %ld\n", 
        // entries[e].cell->id, entries[e].cell->index_entry_pos, entries[e].num_sets);

        COUNT_PARTIAL_INPUT_TIME_END

        entries[e].sets = malloc(sizeof(unsigned long) * entries[e].num_sets);
        entries[e].vector_count = malloc(sizeof(unsigned long) * entries[e].num_sets);
        if (entries[e].sets == NULL)
            exit_with_failure("Error in inv_index.c: Couldn't allocate memory for entry sets.");
        
        COUNT_PARTIAL_INPUT_TIME_START
        fread(entries[e].sets, sizeof(unsigned long), entries[e].num_sets, idx_file);
        fread(entries[e].vector_count, sizeof(unsigned long), entries[e].num_sets, idx_file);
        
        // for(int s = 0; s <  entries[e].num_sets; s++)
        //     printf("set %ld (%ld) - ", entries[e].sets[s], entries[e].vector_count[s]);
        // printf("\n\n");

        COUNT_PARTIAL_INPUT_TIME_END
    }
    

    free(inverted_idx_filename);
    COUNT_PARTIAL_INPUT_TIME_START
    fclose(idx_file);
    COUNT_PARTIAL_INPUT_TIME_END

    printf("(OK)\n");
    return index;
}