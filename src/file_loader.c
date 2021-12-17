#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include "../include/index.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_loader.h"

vector * load_binary_files(unsigned long * curr_total_vectors, const char *bin_files_directory, unsigned int num_files, unsigned int base, unsigned int num_dim)
{
    vector * data_set = NULL;
    unsigned int read_files = 0u; // numbre of read files so far
    unsigned int datasize, table_id, nsets, vector_length;
    uint32_t num_vectors = 0u;
    vector *vector = (struct vector *)malloc(sizeof(struct vector));
    v_type val;
    DIR *dir;
    struct dirent *dfile;
    dir = opendir(bin_files_directory);

    if (dir == NULL)
        exit_with_failure("Error in index.c: Unable to open binary files directory stream!\n");

    vector->values = (v_type *)malloc(sizeof(v_type) * num_dim);

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

            // printf("\nLoading file %s ...\n", dfile->d_name);
            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = "";
            strcat(bin_file_path, bin_files_directory);
            strcat(bin_file_path, "/");
            strcat(bin_file_path, dfile->d_name);

            // get binary table info
            sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize, &table_id, &nsets, &vector_length);

            // check if vector length in file name matches vector length passed as argument
            if (vector_length != num_dim)
                exit_with_failure("Error in file_loader.c:  number of dimentions in index settings does not match vector length in file.\n");

            /* read binary file */
            FILE *bin_file = fopen(bin_file_path, "rb");

            if (bin_file == NULL)
                exit_with_failure("Error in file_loader.c: Binary file not found in directory!\n");

            /* Start processing file: read every vector in binary file */
            int i = 0, j = 0, set_id = 0, total_bytes = base * ((datasize * num_dim) + nsets) / 8;
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
                else if (i <= (unsigned int)num_vectors * num_dim)
                {
                    // end of vector but still in current set
                    if (j > num_dim - 1)
                    {
                        j = 0;
                        // add vector to array of vectors
                        data_set = realloc(data_set, sizeof(struct vector) * (*curr_total_vectors + 1));
                        data_set[*curr_total_vectors].values = malloc(sizeof(v_type) * num_dim);

                        data_set[*curr_total_vectors].set_id = vector->set_id;
                        data_set[*curr_total_vectors].table_id = vector->table_id;
                        vector_cpy(&data_set[*curr_total_vectors], vector, num_dim);

                        *curr_total_vectors = *curr_total_vectors + 1;
                    }

                    fread((void *)(&val), sizeof(val), 1, bin_file);
                    total_bytes -= 4;
                    // printf("\t> value %d: %.10f\n", j, val);
                    vector->values[j] = val;

                    // last value in last vector in current  set
                    if (i == (unsigned int)num_vectors * num_dim)
                    {
                        // add vector to array of vectors
                        data_set = realloc(data_set, sizeof(struct vector) * (*curr_total_vectors + 1));
                        data_set[*curr_total_vectors].values = malloc(sizeof(v_type) * num_dim);

                        data_set[*curr_total_vectors].set_id = vector->set_id;
                        data_set[*curr_total_vectors].table_id = vector->table_id;
                        vector_cpy(&data_set[*curr_total_vectors], vector, num_dim);

                        *curr_total_vectors = *curr_total_vectors + 1;

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
                exit_with_failure("Error in dstree_file_loaders.c: Could not close binary.\n");
        }
    }

    if (read_files == 0)
        exit_with_failure("Error in index.c:  Could not find any binary file in binary files directory.\n");

    free(vector->values);
    free(vector);
    free(dfile);

    return data_set;
}

/* index raw binary vectors (in metric space) */
response index_binary_files(pexeso_index *index, const char *bin_files_directory, unsigned int num_files, unsigned int base)
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
        exit_with_failure("Error in index.c: Unable to open binary files directory stream!\n");

    vector->values = (v_type *)malloc(sizeof(v_type) * index->settings->num_dim);

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
            if (vector_length != index->settings->num_dim)
                exit_with_failure("Error in file_loader.c:  number of dimentions in index settings does not match vector length in file.\n");

            /* read binary file */
            FILE *bin_file = fopen(bin_file_path, "rb");

            if (bin_file == NULL)
                exit_with_failure("Error in file_loader.c: Binary file not found in directory!\n");

            /* Start processing file: read every vector in binary file */
            int i = 0, j = 0, set_id = 0, total_bytes = base * ((datasize * index->settings->num_dim) + nsets) / 8;
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
                else if (i <= (unsigned int)num_vectors * index->settings->num_dim)
                {
                    // end of vector but still in current set
                    if (j > index->settings->num_dim - 1)
                    {
                        j = 0;
                        // insert vector in index
                        if (!insert_vector(index, vector))
                            exit_with_failure("Error in file_loaders.c:  Could not add vector to the index.\n");
                    }

                    fread((void *)(&val), sizeof(val), 1, bin_file);
                    total_bytes -= 4;
                    // printf("\t> value %d: %.10f\n", j, val);
                    vector->values[j] = val;

                    // last value in last vector in current  set
                    if (i == (unsigned int)num_vectors * index->settings->num_dim)
                    {
                        // insert vector in index
                        if (!insert_vector(index, vector))
                            exit_with_failure("Error in file_loader.c:  Could not add vector to the index.\n");

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
                exit_with_failure("Error in dstree_file_loaders.c: Could not close binary.\n");
        }
    }

    if (read_files == 0)
        exit_with_failure("Error in index.c:  Could not find any binary file in binary files directory.\n");

    free(vector->values);
    free(vector);
    free(dfile);

    return OK;
}

bool is_binaryfile(const char *filename) // check if filename has extesion.
{
    char *ext = ".bin";
    size_t nl = strlen(filename), el = strlen(ext);
    return nl >= el && !strcmp(filename + nl - el, ext);
}
