#include <stdio.h>
#include <stdlib.h>
#include "../include/hgrid.h"
#include "../include/level.h"
#include "../include/cell.h"
#include "../include/gsl_matrix.h"
#include "../include/select_pivots.h"
#include "../include/match_map.h"
#include "../include/query_engine.h"
#include "../include/pexeso.h"
#include "../include/file_buffer.h"
#include "../include/file_buffer_manager.h"
#include "../include/file_loader.h"
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include "../include/stats.h"

int main(int argc, char **argv)
{
    /* initialize stats and start measure total time */
    INIT_STATS()
    COUNT_TOTAL_TIME_START

    /* inputs */
    const char * work_dir = "/home/jaouhara/Projects/pexeso-debbug/pexeso/";
    const char * bin_files_directory = "/home/jaouhara/Projects/Dissertation/dssdl/encode/binfiles/"; //target directory
    const char * bin_query_file_directory = "/home/jaouhara/Projects/Dissertation/dssdl/encode/binfiles/query/"; //target directory

    unsigned long total_dl_files = 0ul;

    unsigned int base = 32; // 32 bits to store numbers in binary files
    unsigned int mtr_vector_length = 100, num_dim_metric_space = 100;
    unsigned long long total_vectors = 0ull; // number of vectors in the whole data lake
    unsigned int total_query_vectors = 0u; // query size
    unsigned int max_leaf_size = 76; // max vectors in one leaf cell
    double mtr_buffered_memory_size = 60; // memory  allocated for file buffers (in MB)
    double ps_buffered_memory_size = 10; // memory  allocated for file buffers (in MB)

    double qgrid_mtr_buffered_memory_size = 60; // memory  allocated for file buffers (in MB) for query grid
    double qgrid_ps_buffered_memory_size = 10; // memory  allocated for file buffers (in MB) for query grid

    /* grid settings */
    unsigned int num_levels = 3;  // m
    unsigned int num_pivots = 2;  // number of pivots
    unsigned int fft_scale = 13;   // constant for finding |P| * fft_scale candidate pivots, a good choice of fft_scale is approximately 30 (in paper) and 13 with experiments.
    bool best_fft = 1; // search for the best fft scale that will maximize extremity (have less vectors ou of pivot space)
    unsigned short num_best_fft_iter = 5; // number of iteration to search for best fft
    unsigned int max_fft = 30; // max fft scale (search for the best fft scale that produces the farthest extrimity)
    /* query settings (todo: change threshold to %) */
    float join_threshold = 1; // T 
    float dist_threshold = 0.0; // tau = 0%
    float max_dist = 2; // max euclidean distance between two normalized vectors

    /* num query sets to take from bin_query_file_directory, if equal to -1 the algorithm will take all sets in all files */
    int num_query_sets = -1;
    int min_query_set_size = 0;
    int max_query_set_size = -1;

    unsigned int track_vector = 1; // track vectors id (table_id, column_id)
    unsigned int mode = 0; //  mode 0 = create grid and query
    unsigned int time_to_select_pivots = 0;

    while (1)
    {
        static struct option long_options[] = {
            {"work-dir", required_argument, 0, 'r'},
            {"datalake", required_argument, 0, 'd'},
            {"datalake-size", required_argument, 0, '#'}, // number of files (tables)
            {"queries", required_argument, 0, 'q'},
            {"bits", no_argument, 0, 'b'},

            {"metric-space-dim", required_argument, 0, 'm'},
            {"pivot-space-dim", required_argument, 0, 'p'},

            {"num-levels", required_argument, 0, 'l'},
            {"leaf-size", required_argument, 0, 'i'},
            {"fft-scale", required_argument, 0, 'f'},
            {"best-fft", required_argument, 0, 'e'},
            {"num-best-fft-iter", required_argument, 0, '*'},
            {"max-fft", required_argument, 0, 'k'},
            {"join-threshold", required_argument, 0, 'j'},
            {"dist-threshold", required_argument, 0, 't'},

            {"num-query-sets", required_argument, 0, 'u'},
            {"min-query-set-size", required_argument, 0, 'z'},
            {"max-query-set-size", required_argument, 0, 'a'},

            {"metric-buffer-size", required_argument, 0, 'c'},
            {"pivot-buffer-size", required_argument, 0, 's'},

            {"qgrid-metric-buffer-size", required_argument, 0, '&'},
            {"qgrid-pivot-buffer-size", required_argument, 0, '$'},
            
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

            case '#':
                total_dl_files = atoi(optarg);
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
            
            case 'e':
                best_fft = atoi(optarg);
                if (best_fft != 0 && best_fft != 1)
                {
                    fprintf(stderr, "Please change best fft scale to be 0 or 1.\n");
                    exit(-1);
                }
                break;

            case '*':
                num_best_fft_iter = atoi(optarg);
                break;
            
            case 'k':
                max_fft = atoi(optarg);
                if (max_fft < 0)
                {
                    fprintf(stderr, "Please change max fft scale to be greater than 0.\n");
                    exit(-1);
                }
                break;

            case 'j':
                join_threshold = atof(optarg); // dist threshold is provided in % of the query set size
                if (join_threshold > 1 || join_threshold < 0)
                {
                    fprintf(stderr, "Please change join thershold to be less than 1 and greater than 0.\n");
                    exit(-1);
                }
                break;

            case 't':
                dist_threshold = atof(optarg); // dist threshold is provided in %
                if (dist_threshold > 1 || dist_threshold < 0)
                {
                    fprintf(stderr, "Please change distance thershold to be less than 1 and greater than 0.\n");
                    exit(-1);
                }
                dist_threshold = dist_threshold * max_dist;
                break;

            case 'u':
                num_query_sets = atoi(optarg);
                if (num_query_sets < -1)
                {
                    fprintf(stderr, "Please change number of query sets to be greater than 0.\n");
                    exit(-1);
                }
                break;
            
            case 'z':
                min_query_set_size = atoi(optarg);
                if (min_query_set_size < 0)
                {
                    fprintf(stderr, "Please change min query set size to be greater than or equal to 0.\n");
                    exit(-1);
                }
                break;
            
            case 'a':
                max_query_set_size = atoi(optarg);
                if (max_query_set_size < -1)
                {
                    fprintf(stderr, "Please change max query set size to be greater than or equal to -1.\n");
                    exit(-1);
                }
                break;

            case 'c':
                mtr_buffered_memory_size = atoi(optarg);
                if (mtr_buffered_memory_size < 1)
                {
                    fprintf(stderr, "Please change the bufferd memory size to be greater than 1 MB.\n");
                    exit(-1);
                }
                break;
            
            case '&':
                qgrid_mtr_buffered_memory_size = atoi(optarg);
                if (qgrid_mtr_buffered_memory_size < 1)
                {
                    fprintf(stderr, "Please change the bufferd memory size (for the query grid) to be greater than 1 MB.\n");
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

    /* start count time  to select pivots and make the grid structure */
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

    /* read all vectors in the data set */
    printf("Reading dataset info...");
    if(total_dl_files == 0)
        get_full_datalake_info(bin_files_directory, &total_dl_files, &total_vectors, &mtr_vector_length);
    else
        get_datalake_info(bin_files_directory, total_dl_files, &total_vectors, &mtr_vector_length);
    printf("(OK)\n");

    printf("\tNumber of tables = %lu\n\tNumber of vectors = %llu\n\tVector length in mtric spaces = %u\n\n", total_dl_files, total_vectors, mtr_vector_length);

    printf("Loading dataset files...");
    vector * dataset = load_binary_files(bin_files_directory, 
                                        total_dl_files, total_vectors, base, mtr_vector_length);
    
    if(dataset == NULL)
        exit_with_failure("Error in main.c: Something went wrong, couldn't read dataset vectors!");
    printf("(OK)\n");

    printf("\n\nLooking for pivot vectors...");
    /* search for pivot vectors using pca based algorithm (waiting for response from authors) */
    int dataset_dim [] = {total_vectors, num_dim_metric_space};
    int pivots_mtr_dim [] = {num_pivots, num_dim_metric_space};

    // get pivot vector in from metric dataset
    clock_t start_pivot_selection_time; // start time to select pivots
    clock_t end_pivot_selection_time; // end time after selecting pivots

    vector * pivots_mtr = NULL;
    if(!best_fft) // search for pivots using user defined fft_scale
    {
        COUNT_PIVOT_SELECTION_TIME_START
        pivots_mtr = select_pivots(dataset, dataset_dim, num_pivots, fft_scale);
        COUNT_PIVOT_SELECTION_TIME_END
    }
    else
    {
        COUNT_PIVOT_SELECTION_TIME_START
        pivots_mtr = select_pivots_with_best_fft_scale(dataset, dataset_dim, pivots_mtr_dim, fft_scale, max_fft, num_best_fft_iter);
        COUNT_PIVOT_SELECTION_TIME_END
    }

    // free dataset (dataset was loaded into memory to get pivots)
    for(int dv = total_vectors - 1; dv >=0; dv--)
        free(dataset[dv].values);
    free(dataset);  

    printf("(OK)\n");


    /* map pivot vectors to pivot space pi --> pi' */
    printf("\n\nTransforming pivots to pivot space... ");
    vector * pivots_ps = map_to_pivot_space(pivots_mtr, pivots_mtr_dim, pivots_mtr, num_pivots);
    printf("(OK)\n");
    

    /* pivot space extremity */
    // printf("\nextremity vector (in pivot space):\n");
    vector * pivot_space_extremity = get_extremity(pivots_ps, num_pivots);
    // print_vector(pivot_space_extremity, num_pivots);
    
    /* initialize grid */
    printf("\n\nInitialize grid... ");
    struct query_settings * query_settings = init_query_settings(dist_threshold, join_threshold, num_query_sets, min_query_set_size, max_query_set_size, qgrid_mtr_buffered_memory_size);
    struct grid * grid = (struct grid *) malloc(sizeof(struct grid));
    if (grid == NULL)
        exit_with_failure("Error in main.c: Couldn't allocate memory for grid!");

    if (!init_grid(work_dir, num_pivots, pivots_mtr, pivots_ps, pivot_space_extremity, 
                    num_levels, total_vectors, base, mtr_vector_length, 
                    mtr_buffered_memory_size, max_leaf_size, track_vector, 
                    false, query_settings, grid))
        exit_with_failure("Error in main.c: Couldn't initialize grid!");
    printf("(OK)\n");


    /* initialize grid stats */
    printf("\n\nInitialize grid stats... ");
    if(!init_grid_stats(grid))
        exit_with_failure("Error in main.c: Couldn't initialize grid stats!");
        grid->stats->total_pivot_selection_time += pivot_selection_time;
    printf("(OK)\n");

    /* Build levels */
    printf("\n\nBuilding levels... ");
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

    if (!index_binary_files(grid, index, bin_files_directory, total_dl_files, base))
        exit_with_failure("Error in main.c: Something went wrong, couldn't index binary files.");
    printf("(OK)\n");

    /* count empty leafs in grid */
    int num_empty_leaves = count_empty_leaf_cells(grid);
    grid->stats->empty_leaf_cells_count = num_empty_leaves;

    /* print Rv grid */
    // dump_grid_to_console(grid);

    COUNT_PARTIAL_TIME_END

    grid->stats->grid_building_total_time += partial_time;
    grid->stats->grid_building_input_time += partial_input_time;
    grid->stats->grid_building_output_time += partial_output_time;

    grid->stats->out_of_ps_space_vec_count += out_of_ps_space_vec_count;
    grid->stats->out_of_ps_space_qvec_count += out_of_ps_space_qvec_count;

    grid->stats->loaded_files_count += loaded_files_count;
    grid->stats->loaded_sets_count += loaded_sets_count;
    grid->stats->loaded_files_size += loaded_files_size;
    grid->stats->loaded_vec_count += loaded_vec_count;
    


    /* querying */
    pexeso(bin_query_file_directory, grid, index);


    grid->stats->total_query_time += total_query_time;
    grid->stats->loaded_query_files_count += loaded_query_files_count;
    grid->stats->loaded_query_sets_count += loaded_query_sets_count;
    grid->stats->loaded_query_files_size += loaded_query_files_size;
    grid->stats->loaded_qvec_count += loaded_qvec_count;
    for(int i = 0; i < 7; i++)
        grid->stats->used_lemmas_count[i] = used_lemmas_count[i];


    /* print inverted index */
    // dump_inv_index_to_console(index);

    /* write grid to disk */
    if (!grid_write(grid))
        exit_with_failure("Error main.c:  Could not save the grid to disk.\n");
    grid->stats->grid_building_output_time += partial_output_time;

    /* end of endexing and quering */
    COUNT_TOTAL_TIME_END
    grid->stats->total_time = total_time;

    /* print grid statistics */
    print_grid_stats(grid);


    printf("End of progam: combined indexing and querying time : %.2f secs \n", total_time / 1000000);
    

    /* destroy grid */
    if (!grid_destroy(grid))
        exit_with_failure("Error main.c: Could not destroy grid.\n");
    
    /* destroy inverted index */
    if(!inv_index_destroy(index))
        exit_with_failure("Error main.c: Couldn't destroy inverted index.\n");

    exit(0);
}