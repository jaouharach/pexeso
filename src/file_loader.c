#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include "../include/hgrid.h"
#include "../include/match_map.h"
#include "../include/query_engine.h"

#include "../include/file_loader.h"

vector * load_binary_files(const char *bin_files_directory, unsigned long num_files, unsigned long long total_vectors, unsigned int base, unsigned int mtr_vector_length)
{
    vector * dataset = NULL;
    unsigned int read_files = 0u; // numbre of read files so far
    unsigned long long curr_total_vectors = 0; // numbre of read vectors so far

    // binary file info
    unsigned int datasize, table_id, nsets, vector_length;
    uint32_t num_vectors = 0u;

    v_type val;
    DIR *dir;
    struct dirent *dfile;
    dir = opendir(bin_files_directory);

    if (dir == NULL)
        exit_with_failure("Error in file_loader.c: Unable to open binary files directory stream!\n");

    // temp vector
    vector *vector = (struct vector *)malloc(sizeof(struct vector));
    if (vector == NULL)
        exit_with_failure("Error in file_loader.c: Could not allocate memory for temp vector.");

    vector->values = (v_type *)malloc(sizeof(v_type) * mtr_vector_length);
    if (vector->values == NULL)
        exit_with_failure("Error in file_loader.c: Could not allocate memory for temp vector values.");


    // allocate memory for the whole set of vectors
    dataset = malloc(sizeof(struct vector) * total_vectors);
    if (dataset == NULL)
        exit_with_failure("Error in file_loader.c: Could not allocate memory for dataset.");

    for(int k = 0; k < total_vectors; k++)
    {
        dataset[k].values = malloc(sizeof(v_type) * mtr_vector_length);
        if (dataset[k].values == NULL)
            exit_with_failure("Error in file_loader.c: Could not allocate memory for dataset values.");

    }

    while ((dfile = readdir(dir)) != NULL && num_files > 0) // each file in directory
    {
        if (dfile->d_type != DT_REG) // skip directories
            continue;

        if (is_binaryfile(dfile->d_name))
        {
            num_files--;
            read_files += 1;

            // printf("\nLoading file %s ...\n", dfile->d_name);
            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = "";
            strcat(bin_file_path, bin_files_directory);
            strcat(bin_file_path, "/");
            strcat(bin_file_path, dfile->d_name);

            // get binary table info
            sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize, &table_id, &nsets, &vector_length);

            // check if vector length in file name matches vector length passed as argument
            if (vector_length != mtr_vector_length)
                exit_with_failure("Error in file_loader.c:  number of dimentions in grid settings does not match vector length in file.\n");

            /* read binary file */
            FILE *bin_file = fopen(bin_file_path, "rb");

            if (bin_file == NULL)
                exit_with_failure("Error in file_loader.c: Binary file not found in directory!\n");

            /* Start processing file: read every vector in binary file */
            int i = 0, j = 0, set_id = 0, total_bytes = base * ((datasize * mtr_vector_length) + nsets) / 8;
            // printf("File size in bytes = %u\n\n", total_bytes);

            while (total_bytes)
            {
                if (i == 0)
                {
                    i++;
                    j = 0;
                    //read first integer to check how many vactors in current set
                    fread(&num_vectors, sizeof(num_vectors), 1, bin_file);
                    total_bytes -= 4;
                    // printf("Read set (%u, %u):\n", table_id, set_id);
                    // printf("\t> Num vectors: %u\n", num_vectors);

                    vector->table_id = table_id;
                    vector->set_id = set_id;

                    set_id += 1;
                }
                else if (i <= (unsigned int)num_vectors * mtr_vector_length)
                {
                    // end of vector but still in current set
                    if (j > mtr_vector_length - 1)
                    {
                        j = 0;
                        // add vector to array of vectors
                        dataset[curr_total_vectors].set_id = vector->set_id;
                        dataset[curr_total_vectors].table_id = vector->table_id;
                        dataset[curr_total_vectors].pos = vector->pos;
                        dataset[curr_total_vectors].set_size = vector->set_size;

                        // sanity check: is vector contains nan value exit(1)
                        if(isnan(vector->values[0]))
                            exit_with_failure("Error in file_loader.c: found vector with nan value.");

                        vector_cpy(&dataset[curr_total_vectors], vector, mtr_vector_length);

                        curr_total_vectors = curr_total_vectors + 1;
                    }

                    fread((void *)(&val), sizeof(val), 1, bin_file);
                    total_bytes -= 4;
                    // printf("\t> value %d: %.10f\n", j, val);
                    vector->values[j] = val;

                    // last value in last vector in current  set
                    if (i == (unsigned int)num_vectors * mtr_vector_length)
                    {
                        // add vector to array of vectors
                        dataset[curr_total_vectors].set_id = vector->set_id;
                        dataset[curr_total_vectors].table_id = vector->table_id;
                        dataset[curr_total_vectors].pos = vector->pos;
                        dataset[curr_total_vectors].set_size = vector->set_size;

                        // sanity check: is vector contains nan value exit(1)
                        if(isnan(vector->values[0]))
                            exit_with_failure("Error in file_loader.c: found vector with nan value.");

                        vector_cpy(&dataset[curr_total_vectors], vector, mtr_vector_length);
                        curr_total_vectors = curr_total_vectors + 1;

                        i = 0;
                        j = 0;
                        num_vectors = 0u;
                        continue;
                    }
                    i++;
                    j++;
                }
            }
            if (fclose(bin_file))
                exit_with_failure("Error in file_loader.c: Could not close binary.\n");
        }
    }

    if (read_files == 0)
        exit_with_failure("Error in file_loader.c:  Could not find any binary file in binary files directory.\n");

    // free memory
    free(vector->values);
    free(vector);
    free(dir);

    return dataset;
}

/* index raw binary vectors (in metric space) */
enum response index_binary_files(struct grid *grid, struct inv_index * index, const char *bin_files_directory, unsigned int num_files, unsigned int base)
{
    // printf("step 1.\n");
    unsigned int read_files = 0u; // numbre of read files so far
    unsigned int datasize, table_id, nsets, vector_length;
    uint32_t num_vectors = 0u;
    vector *vector = (struct vector *)malloc(sizeof(struct vector));
    v_type val;
    DIR *dir;
    struct dirent *dfile;
    dir = opendir(bin_files_directory);

    if (dir == NULL)
        exit_with_failure("Error in file_loader.c: Unable to open binary files directory stream!\n");

    vector->values = (v_type *)malloc(sizeof(v_type) * grid->settings->mtr_vector_length);

    if (vector->values == NULL)
        exit_with_failure("Error in file_loader.c: Could not allocate memory for vector values.");

    while ((dfile = readdir(dir)) != NULL && num_files > 0) // each file in directory
    {
        if (dfile->d_type != DT_REG) // skip directories
            continue;

        if (is_binaryfile(dfile->d_name))
        {
            num_files--;
            read_files += 1;

            printf("\n\nIndexing file %s ...\n", dfile->d_name);
            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = "";
            strcat(bin_file_path, bin_files_directory);
            strcat(bin_file_path, "/");
            strcat(bin_file_path, dfile->d_name);

            // get binary table info
            sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize, &table_id, &nsets, &vector_length);

            // check if vector length in file name matches vector length passed as argument
            if (vector_length != grid->settings->mtr_vector_length)
                exit_with_failure("Error in file_loader.c:  number of dimentions in grid settings does not match vector length in file.\n");

            /* read binary file */
            FILE *bin_file = fopen(bin_file_path, "rb");

            if (bin_file == NULL)
                exit_with_failure("Error in file_loader.c: Binary file not found in directory!\n");

            /* Start processing file: read every vector in binary file */
            int i = 0, j = 0, set_id = 0, total_bytes = base * ((datasize * grid->settings->mtr_vector_length) + nsets) / 8;
            printf("File size in bytes = %u\n\n", total_bytes);

            while (total_bytes)
            {
                if (i == 0)
                {
                    i++;
                    j = 0;
                    //read first integer to check how many vactors in current set
                    fread(&num_vectors, sizeof(num_vectors), 1, bin_file);
                    total_bytes -= 4;
                    // printf("Read set (%u, %u):\n", table_id, set_id);
                    // printf("\t> Num vectors: %u\n", num_vectors);

                    vector->table_id = table_id;
                    vector->set_id = set_id;
                    vector->pos = 0; // vector position in set
                    vector->set_size = num_vectors;
                    set_id += 1;
                }
                else if (i <= (unsigned int)num_vectors * grid->settings->mtr_vector_length)
                {
                    // end of vector but still in current set
                    if (j > grid->settings->mtr_vector_length - 1)
                    {
                        j = 0;
                        // printf("vec (%u, %u, %u)\n", vector->table_id, vector->set_id, vector->pos);
                        // insert vector in grid
                        if (!grid_insert(grid, index, vector))
                            exit_with_failure("Error in file_loaders.c:  Could not add vector to the grid.\n");
                        vector->pos = vector->pos + 1;
                    }

                    fread((void *)(&val), sizeof(val), 1, bin_file);
                    total_bytes -= 4;
                    // printf("\t> value %d: %.10f\n", j, val);
                    vector->values[j] = val;

                    // last value in last vector in current  set
                    if (i == (unsigned int)num_vectors * grid->settings->mtr_vector_length)
                    {
                        // printf("vec (%u, %u, %u)\n", vector->table_id, vector->set_id, vector->pos);
                        // insert vector in grid
                        if (!grid_insert(grid, index, vector))
                            exit_with_failure("Error in file_loader.c:  Could not add vector to the grid.\n");

                        i = 0;
                        j = 0;
                        num_vectors = 0u;
                        vector->pos = vector->pos + 1;
                        continue;
                    }
                    i++;
                    j++;
                }
            }
            if (fclose(bin_file))
                exit_with_failure("Error in file_loaders.c: Could not close binary.\n");
        }
    }

    if (read_files == 0)
        exit_with_failure("Error in file_loader.c:  Could not find any binary file in binary files directory.\n");

    // free memory
    free(vector->values);
    free(vector);
    free(dir);

    return OK;
}

/* index raw binary vectors (in metric space), only index a specific number af sets */
struct sid * index_query_binary_files(struct grid *grid, struct grid * Dgrid, struct inv_index * index, const char *bin_files_directory, unsigned int num_files, unsigned int base, int min_query_set_size, int max_query_set_size)
{
    // printf("step 1.\n");
    unsigned int read_files = 0u; // numbre of read files so far
    unsigned int datasize, table_id, nsets, vector_length;

    int num_query_sets = Dgrid->settings->query_settings->num_query_sets;
    unsigned int set_counter = 0;
    struct sid * query_sets = NULL;
    
    uint32_t num_vectors = 0u;
    vector *vector = (struct vector *)malloc(sizeof(struct vector));
    v_type val;
    DIR *dir;
    struct dirent *dfile;
    dir = opendir(bin_files_directory);

    if (dir == NULL)
        exit_with_failure("Error in file_loader.c: Unable to open binary files directory stream!\n");

    vector->values = (v_type *)malloc(sizeof(v_type) * grid->settings->mtr_vector_length);

    if (vector->values == NULL)
        exit_with_failure("Error in file_loader.c: Could not allocate memory for vector values.");

    while ((dfile = readdir(dir)) != NULL && num_files > 0) // each file in directory
    {
        if(num_query_sets == 0)
            break;
        if (dfile->d_type != DT_REG) // skip directories
            continue;

        if (is_binaryfile(dfile->d_name))
        {
            num_files--;
            read_files += 1;

            printf("\n\nIndexing file %s ...\n", dfile->d_name);
            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = "";
            strcat(bin_file_path, bin_files_directory);
            strcat(bin_file_path, "/");
            strcat(bin_file_path, dfile->d_name);

            // get binary table info
            sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize, &table_id, &nsets, &vector_length);

            // check if vector length in file name matches vector length passed as argument
            if (vector_length != grid->settings->mtr_vector_length)
                exit_with_failure("Error in file_loader.c:  number of dimentions in grid settings does not match vector length in file.\n");

            /* read binary file */
            FILE *bin_file = fopen(bin_file_path, "rb");

            if (bin_file == NULL)
                exit_with_failure("Error in file_loader.c: Binary file not found in directory!\n");

            /* Start processing file: read every vector in binary file */
            int i = 0, j = 0, set_id = 0, total_bytes = base * ((datasize * grid->settings->mtr_vector_length) + nsets) / 8;
            printf("File size in bytes = %u\n\n", total_bytes);

            while (total_bytes)
            {
                if (i == 0)
                {
                    if(num_query_sets == 0)
                        break;
                    i++;
                    j = 0;
                    //read first integer to check how many vactors in current set
                    fread(&num_vectors, sizeof(num_vectors), 1, bin_file);
                    total_bytes -= 4;
                    // printf("Read set (%u, %u):\n", table_id, set_id);
                    // printf("\t> Num vectors: %u\n", num_vectors);

                    // skip set if its size doesn't match requirements
                    if(max_query_set_size != -1 || min_query_set_size > 0)
                        if((unsigned int)num_vectors < min_query_set_size || (unsigned int)num_vectors > max_query_set_size)
                        {
                            fseek(bin_file, num_vectors * 4 * grid->settings->mtr_vector_length, SEEK_CUR);
                            i = 0;
                            j = 0;
                            total_bytes -= num_vectors * 4 * vector_length;
                            continue;
                        }
                    // set id for all vectors in the current set
                    vector->table_id = table_id;
                    vector->set_id = set_id;
                    vector->pos = 0; // vector position in set 
                    vector->set_size = num_vectors;
                    
                    // append set id to list of query sets
                    set_counter++;
                    query_sets = realloc(query_sets, sizeof(struct sid) * set_counter);
                    query_sets[set_counter - 1].table_id = table_id;
                    query_sets[set_counter - 1].set_id = set_id;
                    query_sets[set_counter - 1].set_size = num_vectors;

                    if(num_query_sets != -1)
                        num_query_sets--; // read a new set (column)

                    set_id += 1;
                }
                else if (i <= (unsigned int)num_vectors * grid->settings->mtr_vector_length)
                {
                    // end of vector but still in current set
                    if (j > grid->settings->mtr_vector_length - 1)
                    {
                        j = 0;
                        // printf("vec (%u, %u, %u)\n", vector->table_id, vector->set_id, vector->pos);
                        // insert vector in grid
                        if (!grid_insert(grid, index, vector))
                            exit_with_failure("Error in file_loaders.c:  Could not add vector to the grid.\n");
                        vector->pos = vector->pos + 1;
                    }

                    fread((void *)(&val), sizeof(val), 1, bin_file);
                    total_bytes -= 4;
                    // printf("\t> value %d: %.10f\n", j, val);
                    vector->values[j] = val;

                    // last value in last vector in current  set
                    if (i == (unsigned int)num_vectors * grid->settings->mtr_vector_length)
                    {
                        // printf("vec (%u, %u, %u)\n", vector->table_id, vector->set_id, vector->pos);
                        // insert vector in grid
                        if (!grid_insert(grid, index, vector))
                            exit_with_failure("Error in file_loader.c:  Could not add vector to the grid.\n");

                        i = 0;
                        j = 0;
                        num_vectors = 0u;
                        vector->pos = vector->pos + 1;
                        continue;
                    }
                    i++;
                    j++;
                }
            }
            if (fclose(bin_file))
                exit_with_failure("Error in file_loaders.c: Could not close binary.\n");
        }
    }

    if (read_files == 0)
        exit_with_failure("Error in file_loader.c:  Could not find any binary file in binary files directory.\n");

    // free memory
    free(vector->values);
    free(vector);
    free(dir);

    // check if we read more query sets than required
    num_query_sets = Dgrid->settings->query_settings->num_query_sets;
    if(num_query_sets != -1) // if number of query set is specified
        if(set_counter != num_query_sets)
            exit_with_failure("Error in pexeso.c: Function index_query_binary_files() has read more sets than what's required.");
    else
        Dgrid->settings->query_settings->num_query_sets = set_counter;

    return query_sets;
}

bool is_binaryfile(const char *filename) // check if filename has extesion.
{
    char *ext = ".bin";
    size_t nl = strlen(filename), el = strlen(ext);
    return nl >= el && !strcmp(filename + nl - el, ext);
}

/* check dataset directory and count total files and number of vectors */
unsigned long long get_full_datalake_info(const char *bin_files_directory, unsigned long *num_files, unsigned long long *num_vectors, unsigned int *vector_length)
{
    *num_files = 0ul;
    *num_vectors = 0ull;
    *vector_length = 0u;

    // variables found in file name, datasize = total vactors in file, nsets = number of columns/tuples
    unsigned int datasize, table_id, nsets, v_len;

    DIR *dir;
    struct dirent *dfile;
    dir = opendir(bin_files_directory);

    if (dir == NULL)
        exit_with_failure("Error in file_loader.c: Unable to open binary files directory stream!\n");

    // check every file in directory
    while ((dfile = readdir(dir)) != NULL) 
    {
        if (dfile->d_type != DT_REG) // skip directories
            continue;

        if (is_binaryfile(dfile->d_name))
        {
            // get binary file info
            sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize, &table_id, &nsets, &v_len);
            *num_files = *num_files + 1;
            *num_vectors += datasize;

            // check if there exist files with different vector lengths
            if(*vector_length != 0 && *vector_length != v_len)
                exit_with_failure("Error in file_loader.c: Exit dataset directory. Found binary files with different vector lengths!.");

            *vector_length = v_len;
        }
    }
    // free memory
    free(dir);
}

/* check query directory and count total files and number of vectors */
void get_query_data_info(const char *bin_files_directory, int num_query_sets, int min_query_set_size, int max_query_set_size, unsigned long *total_files, unsigned long long *total_vectors, unsigned int *vector_length)
{
    *total_files = 0ul;
    *total_vectors = 0ull;
    *vector_length = 0u;
    int base = 32;

    // variables found in file name, datasize = total vactors in file, nsets = number of columns/tuples
    unsigned int datasize, table_id, nsets, v_len;

    DIR *dir;
    struct dirent *dfile;
    dir = opendir(bin_files_directory);

    uint32_t num_vectors = 0u; // for current set

    if (dir == NULL)
        exit_with_failure("Error in file_loader.c: Unable to open binary files directory stream!\n");

    // check every file in directory
    while ((dfile = readdir(dir)) != NULL) 
    {
        if (dfile->d_type != DT_REG) // skip directories
            continue;

        if (is_binaryfile(dfile->d_name))
        {
            // printf("\n\nReading file %s ...\n", dfile->d_name);
            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = "";
            strcat(bin_file_path, bin_files_directory);
            strcat(bin_file_path, "/");
            strcat(bin_file_path, dfile->d_name);


            // get binary file info
            sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize, &table_id, &nsets, &v_len);
            
            
            /* read binary file */
            FILE *bin_file = fopen(bin_file_path, "rb");
            unsigned int pick_curr_file = 0; // pick current file to be a query file
            if (bin_file == NULL)
                exit_with_failure("Error in file_loader.c: Binary file not found in directory!\n");

            /* Start processing file: read every vector in binary file */
            int i = 0, j = 0, set_id = 0, total_bytes = base * ((datasize * v_len) + nsets) / 8;
            // printf("File size in bytes = %u\n\n", total_bytes);

            while (total_bytes)
            {
                
                if(num_query_sets == 0)
                    break;
                    
                //read first integer to check how many vactors in current set
                fread(&num_vectors, sizeof(num_vectors), 1, bin_file);
                total_bytes -= 4;
                
                // don't count set if its size doesn't match requirements
                if(max_query_set_size != -1 || min_query_set_size > 0)
                    if((unsigned int)num_vectors > min_query_set_size && (unsigned int)num_vectors < max_query_set_size)
                    {
                        *total_vectors += num_vectors;
                        pick_curr_file = 1;
                        if(num_query_sets != -1)
                            num_query_sets--; // read a new set (column)
                    }

                fseek(bin_file, num_vectors * 4 * v_len, SEEK_CUR);
                total_bytes -= num_vectors * 4 * v_len;
            }
            if (fclose(bin_file))
                exit_with_failure("Error in file_loaders.c: Could not close binary.\n");

            if(pick_curr_file == 1)
                *total_files = *total_files + 1;
            *vector_length = v_len;
        }
    }
    // free memory
    free(dir);
}

/* check datalake directory and count number of vectors (number of files is given as input */
unsigned long long get_datalake_info(const char *bin_files_directory, unsigned long num_files, unsigned long long *num_vectors, unsigned int *vector_length)
{
    *num_vectors = 0ull;
    *vector_length = 0u;

    // variables found in file name, datasize = total vactors in file, nsets = number of columns/tuples
    unsigned int datasize, table_id, nsets, v_len;

    DIR *dir;
    struct dirent *dfile;
    dir = opendir(bin_files_directory);

    if (dir == NULL)
        exit_with_failure("Error in file_loader.c: Unable to open binary files directory stream!\n");

    // check every file in directory
    while ((dfile = readdir(dir)) != NULL && num_files > 0) 
    {
        if (dfile->d_type != DT_REG) // skip directories
            continue;

        if (is_binaryfile(dfile->d_name))
        {
            // get binary file info
            sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize, &table_id, &nsets, &v_len);
            num_files = num_files - 1;
            *num_vectors += datasize;

            // check if there exist files with different vector lengths
            if(*vector_length != 0 && *vector_length != v_len)
                exit_with_failure("Error in file_loader.c: Exit dataset directory. Found binary files with different vector lengths!.");

            *vector_length = v_len;
        }
    }
    // free memory
    free(dir);
}

/* save query results to csv file */
enum response  save_to_query_result_file(char * csv_file, struct sid * query_set, struct match_map * map)
{
    if (csv_file == NULL)
        exit_with_failure("Error in file_loader.c: Could not get name of csv file!\n");

    int i,j;
	FILE *fp = NULL;;
	fp = fopen(csv_file,"w+");
    if (fp == NULL)
        exit_with_failure("Error in file_loader.c: Could not open csv file!\n");

	// write header
	fprintf(fp, "TQ:Q, TS:S, qindex, sindex, q, s, d");
	
    // write results
    for(int i = 0; i < map->num_sets; i++){

        if(map->joinable[i] == true)
        {
            fprintf(fp, "\n");
            fprintf(fp,"%u:%u, %u:%u, 0, 0, [], [], na", query_set->table_id, query_set->set_id, map->sets[i].table_id, map->sets[i].set_id);
        }
    }
    fclose(fp);

    return OK;
}

/* make result file name and path */
char * make_file_path(char * work_dir, struct sid * query_set, unsigned int l, unsigned int dlsize, unsigned int vector_length, float runtime, unsigned int total_checked_vec)
{
    DIR* dir = opendir(work_dir);
	if (!dir)
    {
		printf("WARNING! Experiment direstory '%s' does not exist!", work_dir);
		exit(1);
	}
    int string_size = get_ndigits(query_set->table_id) + get_ndigits(query_set->set_id) + get_ndigits(l)
							 + get_ndigits(dlsize) + get_ndigits(vector_length) + get_ndigits((unsigned int) runtime) + get_ndigits(total_checked_vec)
							 + get_ndigits(query_set->set_size) + strlen("TQ_Q_qsize_l_dlsize_len_runtime_ndistcalc_dataaccess.csv")
							 + strlen(work_dir)
							 + 4 // float decimal precision for dlsize and runtime (.00)
							 + 5;

    char * filepath = malloc(string_size);

	sprintf(filepath, "%s/TQ%u_Q%u_qsize%u_l%u_dlsize%u_len%u_runtime%.2f_ndistcalc_dataaccess%u.csv"
			, work_dir, query_set->table_id, query_set->set_id, query_set->set_size, l, dlsize, vector_length, runtime, total_checked_vec);
    
    filepath[string_size - 1] = '\0';
    
    free(dir);
	return filepath;
}

/* create directory to store query results */
char * make_result_directory(char * work_dir, char* algorithm, unsigned int l, unsigned int num_query_sets, int min_query_set_size, int max_query_set_size)
{
	char * result_dir_name = malloc(get_ndigits(l) + get_ndigits(num_query_sets)
									+ get_ndigits(min_query_set_size) + get_ndigits(max_query_set_size)
									+ strlen("/_l_q_min_max") + strlen(work_dir)+ strlen(algorithm) + 1);

	sprintf(result_dir_name, "%s/%s_l%u_%uq_min%d_max%d", work_dir, algorithm, l, num_query_sets, min_query_set_size, max_query_set_size);

	printf("Result directory name: %s\n", result_dir_name);
	DIR* dir = opendir(result_dir_name);
	if (dir)
  {
      printf("WARNING! Results directory already exists. Please delete directory : %s.\n", result_dir_name);
      exit(-1);
  }
  mkdir(result_dir_name, 0777);
  
  return result_dir_name;
}

/* save query results to disk */
enum response save_results_to_disk(struct grid * Dgrid, struct grid * Qgrid, struct match_map * match_map)
{
    char * algorithm = "pexeso";
    char * work_dir = Dgrid->settings->work_directory;
    unsigned int dataset_size = Dgrid->total_records;
    unsigned int mtr_vector_length = Dgrid->settings->mtr_vector_length;
    unsigned int l = 0; // total tables indexed in grid
    unsigned int dlsize = 0; // dataset size in MB
    int num_query_sets = Dgrid->settings->query_settings->num_query_sets;
    int min_query_set_size = Dgrid->settings->query_settings->min_query_set_size;
    int max_query_set_size = Dgrid->settings->query_settings->max_query_set_size;

    char * result_dir = make_result_directory(work_dir, algorithm, l, num_query_sets, min_query_set_size, max_query_set_size);
    if(result_dir == NULL)
        exit_with_failure("Error in file_loader.c: Couldn't make result directory.");

    // loop match maps for all query sets
    struct sid *query_set = NULL;
    struct match_map *curr_map;
    unsigned int total_checked_vectors = 0;
    float runtime = 0.0;

    for(int m = 0; m < num_query_sets; m++)
    {
        curr_map = &match_map[m];
        query_set = &curr_map->query_set;
        total_checked_vectors = curr_map->total_checked_vectors;
        runtime = curr_map->query_time;

        char * file_path = make_file_path(result_dir, query_set, l, dlsize, mtr_vector_length, runtime, total_checked_vectors);
        
        if(!save_to_query_result_file(file_path, query_set, curr_map))
            exit_with_failure("Error in file_loader.c: Couldn't save query results to csv file.");
        
        free(file_path);
    }
    
    free(result_dir);
    return OK;
}