#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/file_loader.h"
#include "../include/gsl_matrix.h"
#include "../include/select_pivots.h"
#include "../include/file_buffer.h"
#include "../include/file_buffer_manager.h"
#include "../include/match_map.h"
#include "../include/query_engine.h"
#include "../include/pexeso.h"
#include <unistd.h>
#include <getopt.h>


int main(int argc, char **argv)
{
    /* dataset */
    const char * work_dir = "/home/jaouhara/Projects/pexeso-debbug/pexeso/Hgrid/";
    const char * bin_files_directory = "/home/jaouhara/Projects/Dissertation/dssdl/encode/binfiles/"; //target directory
    const char * bin_query_file_directory = "/home/jaouhara/Projects/Dissertation/dssdl/encode/binfiles/query/"; //target directory

    unsigned long num_files = 0ul;
    unsigned int base = 32; // 32 bits to store numbers in binary files
    unsigned int mtr_vector_length = 100, num_dim_metric_space = 100;
    unsigned long long total_vectors = 0ull; // number of vectors in the whole data lake
    unsigned int total_query_vectors = 0u; // query size
    unsigned int max_leaf_size = 76; // max vectors in one leaf cell
    double mtr_buffered_memory_size = 60; // memory  allocated for file buffers (in MB)
    double ps_buffered_memory_size = 10; // memory  allocated for file buffers (in MB)

    /* grid settings */
    unsigned int num_levels = 3;  // m
    unsigned int num_pivots = 2;  // number of pivots
    unsigned int fft_scale = 13;   // constant for finding |P| * fft_scale candidate pivots, a good choice of fft_scale is approximately 30 (in paper) and 13 with experiments.

    /* query settings (todo: change threshold to %) */
    unsigned int join_threshold = 5; // T 
    float dist_threshold = 0.25; // tau


    unsigned int track_vector = 1; // track vectors id (table_id, column_id)
    unsigned int mode = 0; //  mode 0 = create grid and query

    while (1)
    {
        static struct option long_options[] = {
            {"work-dir", required_argument, 0, 'r'},
            {"dataset", required_argument, 0, 'd'},
            {"queries", required_argument, 0, 'q'},
            {"bits", no_argument, 0, 'b'},

            {"metric-space-dim", required_argument, 0, 'm'},
            {"pivot-space-dim", required_argument, 0, 'p'},

            {"num-levels", required_argument, 0, 'l'},
            {"leaf-size", required_argument, 0, 'i'},
            {"fft-scale", required_argument, 0, 'f'},
            {"join-threshold", required_argument, 0, 'j'},
            {"dist-threshold", required_argument, 0, 't'},
            
            {"metric-buffer-size", required_argument, 0, 'c'},
            {"pivot-buffer-size", required_argument, 0, 's'},
            
            {"track-vector", no_argument, 0, 'v'},
            {"mode", required_argument, 0, 'x'},
            {"help", no_argument, 0, '?'}
        };

        /* store the option index. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "", long_options, &option_index);

        if (c == -1)
            break;
        switch (c)
        {
            case 'r':
                work_dir = optarg;
                break;

            case 'd':
                bin_files_directory = optarg;
                break;

            case 'q':
                bin_query_file_directory = optarg;
                break;

            case 'b':
                base = atof(optarg);
                break;

            case 'm':
                mtr_vector_length = atoi(optarg);
                num_dim_metric_space = atoi(optarg);
                if (mtr_vector_length < 0)
                {
                    fprintf(stderr, "Please change metric space dimension to be greater than 0.\n");
                    exit(-1);
                }
                break;

            case 'p':
                num_pivots = atoi(optarg);
                if (num_pivots < 0)
                {
                    fprintf(stderr, "Please change pivot space dimension to be greater than 0.\n");
                    exit(-1);
                }
                break;

            case 'l':
                num_levels = atoi(optarg);
                if (num_levels < 0)
                {
                    fprintf(stderr, "Please change number of levels to be greater than 0.\n");
                    exit(-1);
                }
                break;
            
            case 'i':
                max_leaf_size = atoi(optarg);
                if (max_leaf_size < 0)
                {
                    fprintf(stderr, "Please change leaf size to be greater than 0.\n");
                    exit(-1);
                }
                break;

            case 'f':
                fft_scale = atoi(optarg);
                if (fft_scale < 0)
                {
                    fprintf(stderr, "Please change fft scale to be greater than 0.\n");
                    exit(-1);
                }
                break;
            
            case 'j':
                join_threshold = atoi(optarg);
                //  (todo) change join threshold to be in %
                // if (join_threshold > 1 || join_threshold < 0)
                // {
                //     fprintf(stderr, "Please change join thershold to be less than 1 and greater than 0.\n");
                //     exit(-1);
                // }
                break;

            case 't':
                dist_threshold = atof(optarg);
                //  (todo) change distance threshold to be in %
                // if (dist_threshold > 1 || dist_threshold < 0)
                // {
                //     fprintf(stderr, "Please change distance thershold to be less than 1 and greater than 0.\n");
                //     exit(-1);
                // }
                break;

            case 'c':
                mtr_buffered_memory_size = atoi(optarg);
                if (mtr_buffered_memory_size < 10)
                {
                    fprintf(stderr, "Please change the bufferd memory size to be greater than 10 MB.\n");
                    exit(-1);
                }
                break;

            case 's':
                ps_buffered_memory_size = atoi(optarg);
                if (ps_buffered_memory_size < 10)
                {
                    fprintf(stderr, "Please change the bufferd memory size to be greater than 10 MB.\n");
                    exit(-1);
                }
                break;

            case 'v':
                track_vector = atoi(optarg);
                break;

            case 'x':
                mode = atoi(optarg);
                break;

            case '?':
                printf("Usage:\n\
                            \t--dataset XX \t\t\tThe path to the dataset directory\n\
                            \t--queries XX \t\t\tThe path to the queries directory\n\
                            \t--dataset-size XX \t\tThe number of time series to load\n\
                            \t--queries-size XX \t\tThe number of queries to run\n\
                            \t--mode: 0 = index & query (only option)\t\t\n");

                return 0;
                break;

            default:
                exit(-1);
                break;
        }
    }

    /* read all vectors in the data set */
    printf("Reading dataset info...");
    get_dataset_info(bin_files_directory, &num_files, &total_vectors, &mtr_vector_length);
    printf("(OK)\n");

    printf("\tNumber of (.bin) files = %lu\n\tNumber of vectors = %llu\n\tVector length in mtric spaces = %u\n\n", num_files, total_vectors, mtr_vector_length);

    printf("Loading dataset files...");
    vector * dataset = load_binary_files(bin_files_directory, 
                                        num_files, total_vectors, base, mtr_vector_length);
    
    // printf("Dataset\n");
    // for(int v = 0; v < total_vectors; v++)
    // {
    //     printf("\nvid:(%d, %d)", dataset[v].table_id, dataset[v].set_id);
    //     print_vector(&dataset[v], mtr_vector_length);
    // }
    // printf("End of dataset\n");


    if(dataset == NULL)
        exit_with_failure("Error in main.c: Something went wrong, couldn't read dataset vectors!");
    printf("(OK)\n");

    printf("\n\nLooking for pivot vectors...");
    /* search for pivot vectors using pca based algorithm (waiting for response from authors) */
    int dataset_dim [] = {total_vectors, num_dim_metric_space};
    int pivots_mtr_dim [] = {num_pivots, num_dim_metric_space};

    // get pivot vector in from metric dataset
    vector * pivots_mtr = select_pivots(dataset, dataset_dim, num_pivots, fft_scale);
    
    
    // printf("selected pivots (metric space)\n");
    // for(int p = 0; p < num_pivots; p++)
    //     print_vector(&pivots_mtr[p], num_dim_metric_space);


    // free dataset (dataset was loaded into memory to get pivots)
    for(int dv = total_vectors - 1; dv >=0; dv--)
        free(dataset[dv].values);
    free(dataset);  

    printf("(OK)\n");


    /* map pivot vectors to pivot space pi --> pi' */
    printf("\n\nTransforming pivots to pivot space... ");
    vector * pivots_ps = map_to_pivot_space(pivots_mtr, pivots_mtr_dim, pivots_mtr, num_pivots);

    // printf("selected pivots (pivot space)\n");
    // for(int p = 0; p < num_pivots; p++)
    //     print_vector(&pivots_ps[p], num_pivots);

    printf("(OK)\n");
    

    /* pivot space extremity */
    printf("\nextremity vector (in pivot space):\n");
    // vector * pivot_space_extremity = get_rand_vector(num_pivots);
    vector * pivot_space_extremity = get_extremity(pivots_ps, num_pivots);
    print_vector(pivot_space_extremity, num_pivots);
    


    /* initialize grid */
    printf("\n\nInitialize grid... ");
    struct query_settings * query_settings = init_query_settings(dist_threshold, join_threshold);

    struct grid * grid = (struct grid *) malloc(sizeof(struct grid));
    if (grid == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for grid!");

    if (!init_grid(work_dir, num_pivots, pivots_mtr, pivots_ps, pivot_space_extremity, 
                    num_levels, total_vectors, base, mtr_vector_length, 
                    mtr_buffered_memory_size, ps_buffered_memory_size, max_leaf_size, track_vector, 
                    false, query_settings, grid))
        exit_with_failure("Error in main.c: Couldn't initialize grid!");

    printf("(OK)\n");

    

    /* Display settings */
    printf("\n\t\t*** \tGRID SETTINGS\t ***\t\n");
    printf("--------------------------------------------------------------\n");
    printf("\t\tNumber of pivots = %d\n", grid->settings->num_pivots);
    printf("\t\tNumber of levels = %d\n", grid->settings->num_levels);
    printf("\t\tPivot space volume = %f\n", grid->settings->pivot_space_volume);
    printf("\t\tNumber of leaf cells = %d\n", grid->settings->num_leaf_cells);
    printf("\t\tLeaf cell edge length = %f\n", grid->settings->leaf_cell_edge_length);
    
    

    printf("\t\tPivot vectors (in pivot space):\n\n");
    for(int i = 0; i < num_pivots; i++)
    {
        print_vector(&grid->settings->pivots_ps[i], num_pivots);
    }
    printf("--------------------------------------------------------------\n");


    /* Build levels */
    printf("\n\nBuild levels... ");
    if(!init_root(grid))
        exit_with_failure("Error in main.c: Couldn't initialize root level!");

    if (!init_first_level(grid))
        exit_with_failure("Error in main.c: Couldn't initialize first level!");

    if (!init_levels(grid))
        exit_with_failure("Error in main.c: Couldn't initialize grid levels!");
    printf("(OK)\n");

    
    /* insert dataset in grid (read and index all vectors in the data set) */
    printf("Index dataset vectors and build inverted index...");
    struct inv_index * index = malloc(sizeof(struct inv_index));
    index->num_entries = 0; index->num_distinct_sets = 0;
    
    if(index == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for inverted index.");

    if (!index_binary_files(grid, index, bin_files_directory, num_files, base))
        exit_with_failure("Error in main.c: Something went wrong, couldn't index binary files.");
    printf("(OK)\n");

    printf("\n\nBuild match and mismatch map...");
    struct match_map * match_map = malloc(sizeof(struct match_map));
    if(match_map == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for match map.");

    if (!init_match_map(index, match_map))
        exit_with_failure("Error in main.c: Couldn't initialize match map!");
    printf("(OK)\n");

    /* print Rv grid */
    dump_grid_to_console(grid);

    /* querying */
    pexeso(bin_query_file_directory, grid, index, match_map);

    /* print inverted index */
    // dump_inv_index_to_console(index);

    /* print match map */
    dump_match_map_to_console(match_map);

    /* write grid to disk */
    if (!grid_write(grid))
        exit_with_failure("Error main.c:  Could not save the grid to disk.\n");
    
    /* destroy grid */
    if (!grid_destroy(grid))
        exit_with_failure("Error main.c: Could not destroy grid.\n");
    
    /* destroy inverted index */
    if(!inv_index_destroy(index))
        exit_with_failure("Error main.c: Couldn't destroy inverted index.\n");

    /* destroy match map */
    if(!match_map_destroy(match_map))
        exit_with_failure("Error main.c: Couldn't destroy match map.\n");
    
    exit(0);
}