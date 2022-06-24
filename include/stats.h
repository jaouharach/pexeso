// benchmark
#include <time.h>
#include <sys/time.h>

// counters
unsigned long total_cells_count;
unsigned long leaf_cells_count;
unsigned long empty_leaf_cells_count;

unsigned long loaded_vec_count;
unsigned long out_of_ps_space_vec_count;
unsigned long out_of_ps_space_qvec_count;

unsigned long checked_cells_count;

unsigned long total_queries_count;

// timers 
double start;
double end;

struct timeval total_time_start; // query time and index creation
struct timeval current_time;


struct timeval partial_time_start; // time to select pivots, create the grid and index files
struct timeval partial_input_time_start;
struct timeval partial_output_time_start;
struct timeval pivot_selection_time_start;

double total_input_time;
double total_output_time;
double total_parse_time;
double total_query_time;

double total_time;

double partial_time;
double partial_input_time;
double partial_output_time;

double pivot_selection_time;



#define INIT_STATS() total_time = 0;\
                    total_input_time = 0;\
                    total_output_time = 0;\
                    total_query_time = 0;\
                    total_parse_time = 0;\
                    pivot_selection_time = 0;\
                    total_cells_count = 0;\
                    leaf_cells_count = 0;\
                    empty_leaf_cells_count =0;\
                    loaded_vec_count = 0;\
                    checked_cells_count = 0;\
                    out_of_ps_space_vec_count = 0;\
                    out_of_ps_space_qvec_count = 0;

#define COUNT_NEW_OUT_OF_PIVOT_SPACE_VECTOR ++out_of_ps_space_vec_count;
#define COUNT_NEW_OUT_OF_PIVOT_SPACE_QUERY_VECTOR ++out_of_ps_space_qvec_count;

#define RESET_PARTIAL_COUNTERS() partial_time = 0;\
                                partial_input_time = 0;\
				                partial_output_time = 0;\

#define COUNT_TOTAL_TIME_START gettimeofday(&total_time_start, NULL);
#define COUNT_TOTAL_TIME_END  gettimeofday(&current_time, NULL); \
                                start = total_time_start.tv_sec * 1000000 + (total_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000  + (current_time.tv_usec); \
                                total_time += (end - start);

#define COUNT_PIVOT_SELECTION_TIME_START gettimeofday(&pivot_selection_time_start, NULL);

#define COUNT_PIVOT_SELECTION_TIME_END  gettimeofday(&current_time, NULL); \
                                start = pivot_selection_time_start.tv_sec * 1000000 + (pivot_selection_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000  + (current_time.tv_usec); \
                                pivot_selection_time += (end - start);

#define COUNT_PARTIAL_TIME_START gettimeofday(&partial_time_start, NULL);
#define COUNT_PARTIAL_INPUT_TIME_START gettimeofday(&partial_input_time_start, NULL);   
#define COUNT_PARTIAL_OUTPUT_TIME_START gettimeofday(&partial_output_time_start, NULL);
        
#define COUNT_PARTIAL_TIME_END  gettimeofday(&current_time, NULL); \
                                start = partial_time_start.tv_sec * 1000000 + (partial_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000  + (current_time.tv_usec); \
                                partial_time += (end - start);

#define COUNT_PARTIAL_INPUT_TIME_END  gettimeofday(&current_time, NULL); \
                                start = partial_input_time_start.tv_sec * 1000000 + (partial_input_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000 + (current_time.tv_usec); \
                                partial_input_time += (end - start);

#define COUNT_PARTIAL_OUTPUT_TIME_END gettimeofday(&current_time, NULL); \
                                start = partial_output_time_start.tv_sec * 1000000 + (partial_output_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000  + (current_time.tv_usec); \
                                partial_output_time += (end - start);
