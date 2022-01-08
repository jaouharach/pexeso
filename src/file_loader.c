#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "../include/hgrid.h"
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
response index_binary_files(struct grid *grid, struct inv_index * index, const char *bin_files_directory, unsigned int num_files, unsigned int base)
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

                    set_id += 1;
                }
                else if (i <= (unsigned int)num_vectors * grid->settings->mtr_vector_length)
                {
                    // end of vector but still in current set
                    if (j > grid->settings->mtr_vector_length - 1)
                    {
                        j = 0;
                        // insert vector in grid
                        if (!grid_insert(grid, index, vector))
                            exit_with_failure("Error in file_loaders.c:  Could not add vector to the grid.\n");
                    }

                    fread((void *)(&val), sizeof(val), 1, bin_file);
                    total_bytes -= 4;
                    // printf("\t> value %d: %.10f\n", j, val);
                    vector->values[j] = val;

                    // last value in last vector in current  set
                    if (i == (unsigned int)num_vectors * grid->settings->mtr_vector_length)
                    {
                        // insert vector in grid
                        if (!grid_insert(grid, index, vector))
                            exit_with_failure("Error in file_loader.c:  Could not add vector to the grid.\n");

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

bool is_binaryfile(const char *filename) // check if filename has extesion.
{
    char *ext = ".bin";
    size_t nl = strlen(filename), el = strlen(ext);
    return nl >= el && !strcmp(filename + nl - el, ext);
}

/* check dataset directory and count total files and number of vectors */
unsigned long long get_dataset_info(const char *bin_files_directory, unsigned long *num_files, unsigned long long *num_vectors, unsigned int *vector_length)
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