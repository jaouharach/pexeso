#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct vector vector;

struct vector
{
    float *values;
};

int main(int argc, char const *argv[])
{
    int min = 1, max = 5, num_dim = 2, i, j;
    float len = 1.0; //cell length

    // check if num cells is integer
    if (fmod(pow(max - min, num_dim), pow(len, num_dim)) == 0)
    {
        int nc = pow(max - min, num_dim) / pow(len, num_dim); //number of cells
        int len_row = sqrt(nc);                               // number of cells in one row

        float *coor_values = (float *)malloc(sizeof(float) * len_row); // list of all possible (distinct) coordinate values

        printf("Number of cells = %d\t row length = %d\n", nc, len_row);

        for (i = 0, j = 1; i < len_row; i++, j += 2)
        {
            coor_values[i] = max - ((j / 2) * len);
        }

        // init centers => allocate memory
        vector *centers = (struct vector *)malloc(sizeof(struct vector) * nc);
        for (i = 0; i < nc; i++)
        {
            centers[i].values = (float *)malloc(sizeof(float) * num_dim);
        }

        /* centers is a list of all possible combination of coor_values */
        
    }
    else
        printf("Warning: Number of cells can only be integer! please change settings.");

    return 0;
}

// printf("Dinstinct coordinate values.\n");
// for(i = 0; i < len_row; i++)
// {
//     printf("v%d = %.2f\n", i+1, coor_values[i]);
// }