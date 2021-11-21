#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct vector
{
    int table_id;
    int set_id;
    float *values;
} vector;

void find_combinations(int V[], int d, int dim, vector *vectors, vector out, int append_at);

static int curr_vector = 0;

int main()
{
    int d = 3;              // N dimentions
    int V[] = {1, 2, 3, 4}; // distinct coordinate values for center points
    int len_v = 4;
    int n_vectors = pow(len_v, d); // number of possible combinations
    int p, x;

    printf("num vectors = %d\n", n_vectors);

    // init list of vectors
    vector out;
    vector *vectors = malloc(sizeof(struct vector) * n_vectors);

    out.values = malloc(sizeof(float) * d);
    for (p = 0; p < n_vectors; p++)
    {
        vectors[p].values = malloc(sizeof(v_type) * d);
    }

    find_combinations(V, d, d, vectors, out, 0);

    for(p = 0; p < n_vectors; p++)
    {
        printf("Vector %d:\n(", p+1);
        for(x = 0; x < d; x++)
        {
            printf("%.2f, ", vectors[p].values[x]);
        }
        printf(")\n\n");
    }


    free(vectors);
    free(out.values);

    return 0;
}

void find_combinations(int V[], int d, int dim, vector * vectors, vector out, int append_at)
{
    int len_v = 4;
    if (len_v == 0 || d > len_v)
    {
        printf("Error!\n");
        return;
    }

    if (d == 0)
    {
        // add out to sub_array
        printf("d = 0.\n\n");
        for(int x = 0; x < dim; x++)
            vectors[curr_vector].values[x] = out.values[x];
        curr_vector++;
        return;
    }

    int j;
    for (j = 0; j < len_v; j++)
    {
        // add V[j] to out values
        out.values[append_at] = V[j];
        find_combinations(V, d - 1, dim, vectors, out, append_at + 1);;
    }

    
}