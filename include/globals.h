#include<stdio.h>

#define OK 1;
#define FAILED 0;

#ifndef GUARD_DES_Utils /* prevents errors when including twice */
#define GUARD_DES_Utils
#endif
typedef struct index index;
typedef struct index_settings index_settings;
typedef struct level level;
typedef struct cell cell;
typedef struct vector vector;

void exit_with_error(char * message);
