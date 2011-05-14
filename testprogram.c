#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "sse3d.h"

void print_matrix(const char *name, const sse3d_matrix_t *matrix)
{
    int i;

    printf("MATRIX: %s\n", name);
    for (i = 0; i < 16; i++)
    {
        if (i%4 == 0) printf("[");
        printf("%0.3f ", matrix->_v[i]);
        if (i%4 == 3) printf("]\n");
    }
    printf("\n");
}

void print_vector(const char *name, const sse3d_vector_t *vector)
{
    printf("VECTOR: %s\n[%0.3f %0.3f %0.3f %0.3f]\n\n", name, vector->x, vector->y, vector->z, vector->w);
}

int main()
{
    aligned sse3d_matrix_t matrix;
    aligned sse3d_matrix_t translation;
    aligned sse3d_matrix_t rotation;
    aligned sse3d_matrix_t scale;
    aligned sse3d_matrix_t lookat;

    aligned sse3d_vector_t vector;
    aligned sse3d_vector_t target   = { 0.0f, -4.0f, -3.0f, 0.0f };
    aligned sse3d_vector_t camera   = { 0.0f,  0.0f,  1.0f, 0.0f };
    aligned sse3d_vector_t up_vec   = { 0.0f,  1.0f,  0.0f, 0.0f };
    aligned sse3d_vector_t original = { 1.0f,  2.0f,  3.0f, 1.0f };

    int i;
    int nr_matrix_tests = 10000000;
    int nr_vec = 30000000;
    sse3d_vector_t* vector_list;
    clock_t start, end;
    

    sse3d_translation_matrix(&translation, 3, 3, 0);
    sse3d_rotation_z_matrix(&rotation, M_PI/4);
    sse3d_scale_matrix(&scale, 10, 10, 10);
    sse3d_lookat_matrix(&lookat, &camera, &target, &up_vec);

    print_matrix("translation", &translation);
    print_matrix("rotation", &rotation);
    print_matrix("scale", &scale);
    print_matrix("lookat", &lookat);

    memcpy(&matrix, &rotation, sizeof(sse3d_matrix_t));

    sse3d_multiply_matrix(&matrix, &matrix, &matrix);
    sse3d_multiply_matrix(&matrix, &matrix, &matrix);
    sse3d_multiply_matrix(&matrix, &matrix, &matrix);

    print_matrix("rotation^8 (identity)", &matrix);
    
    print_vector("original", &original);
    sse3d_multiply_vectors(&vector, &translation, &original, 1);
    print_vector("translated", &vector);
    sse3d_multiply_vectors(&vector, &rotation, &vector, 1);
    print_vector("translated * rotated", &vector);
    sse3d_multiply_vectors(&vector, &scale, &vector, 1);
    print_vector("translated * rotated * scaled", &vector);
    
    sse3d_multiply_matrix(&matrix, &scale, &rotation);
    sse3d_multiply_matrix(&matrix, &matrix, &translation);
    sse3d_multiply_vectors(&vector, &matrix, &original, 1);
    print_vector("translated * rotated * scaled (via matrix premultiply)", &vector);

    vector_list = aligned_malloc(sizeof(sse3d_vector_t) * nr_vec);
    
    start = clock();
    sse3d_multiply_vectors(vector_list, &lookat, vector_list, nr_vec);
    end = clock();

    printf("Vertices per second: %0.0f\n", 1.0f / ((end-start) / (float)CLOCKS_PER_SEC) * nr_vec);

    aligned_free(vector_list);

    start = clock();
    for (i=0; i<nr_matrix_tests; i++)
    {
        sse3d_multiply_matrix(&matrix, &matrix, &matrix);
    }
    end = clock();

    printf("Matrix muliplications per second: %0.0f\n", 1.0f / ((end-start) / (float)CLOCKS_PER_SEC) * nr_matrix_tests);

    getchar();
    return 0;
}