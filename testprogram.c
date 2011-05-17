#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "sse3d.h"
#include "sse3d_window.h"

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

aligned sse3d_vector_t vertices[10000];
aligned sse3d_vector_t normals[10000];
aligned sse3d_vector_t v_transformed[10000];
aligned sse3d_vector_t n_transformed[10000];
int indices[30000];

int nr_vertices, nr_indices, nr_normals;

void render(unsigned short *z_buffer, unsigned int *n_buffer, int width, int height)
{
    int i;
    static float angle = 0;
    
    static aligned sse3d_matrix_t model, model_scale, model_rotation, model_rotation_y, model_rotation_x,model_translation;
    static aligned sse3d_matrix_t projection, projection_scale, projection_translation, lookat;
    static aligned sse3d_matrix_t identity, transform;

    sse3d_identity_matrix(&identity);
    sse3d_scale_matrix(&model_scale, .02f, .02f, .02f);
    sse3d_rotation_x_matrix(&model_rotation, -M_PI/2.0f);
    sse3d_rotation_x_matrix(&model_rotation_x, 0);
    sse3d_rotation_y_matrix(&model_rotation_y, angle += 0.01f);
    sse3d_translation_matrix(&model_translation, 0.0f, -0.1f, 0.15f);
    
    sse3d_multiply_matrix(&model, &model_rotation_y, &model_rotation);
    sse3d_multiply_matrix(&model, &model_rotation_x, &model);
    sse3d_multiply_matrix(&model, &model_scale, &model);
    sse3d_multiply_matrix(&model, &model_translation, &model);
    
    sse3d_scale_matrix(&projection_scale, 2.7f * height / 2.0f, 2.7f * height / 2.0f, 1.0f);
    sse3d_translation_matrix(&projection_translation, width / 2.0f, height / 2.0f, 0.0f);
    sse3d_multiply_matrix(&projection, &projection_translation, &projection_scale);

    sse3d_multiply_matrix(&transform, &projection, &model);
    sse3d_multiply_vectors(v_transformed, &transform, vertices, nr_vertices);
    sse3d_multiply_vectors(n_transformed, &model, normals, nr_normals);

    sse3d_prepare_render_vectors(v_transformed, nr_vertices);
    for (i=0; i<nr_indices; i+=6)
    {
        sse3d_draw_triangle(z_buffer, n_buffer, width, height, 
            v_transformed + indices[i+0] - 1, v_transformed + indices[i+1] - 1, v_transformed + indices[i+2] - 1,
            n_transformed + indices[i+3] - 1, n_transformed + indices[i+4] - 1, n_transformed + indices[i+5] - 1);
    }
}


int WINAPI WinMain(HINSTANCE instance, HINSTANCE prevInstance, LPSTR cmdLine, int showCmd)
{
    MSG msg;
    RECT windowRect;
    char buffer[1024];

    sse3d_window_t* window;

    FILE *dragon = fopen("dragon.txt", "r");
    while (fgets(buffer, 1023, dragon))
    {
        if (buffer[0] == 'v' && buffer[1] == ' ')
        {
            sscanf(buffer, "v %f %f %f", &vertices[nr_vertices].x, &vertices[nr_vertices].y, &vertices[nr_vertices].z);
            vertices[nr_vertices].w = 1;
            nr_vertices++;
        }
        if (buffer[0] == 'v' && buffer[1] == 'n')
        {
            sscanf(buffer, "vn %f %f %f", &normals[nr_normals].x, &normals[nr_normals].y, &normals[nr_normals].z);
            normals[nr_normals].w = 0;
            nr_normals++;
        }
        if (buffer[0] == 'f' && buffer[1] == ' ')
        {
            int a, b, c, na, nb, nc;
            if (sscanf(buffer, "f %d/%*d/%d %d/%*d/%d %d/%*d/%d", &a, &na, &b, &nb, &c, &nc))
            {
                indices[nr_indices++] = a;
                indices[nr_indices++] = b;
                indices[nr_indices++] = c;
                indices[nr_indices++] = na;
                indices[nr_indices++] = nb;
                indices[nr_indices++] = nc;
            }
        }
    }

    window = sse3d_create_window(instance, render);
    SetRect(&windowRect, 0, 0, FRAMEBUFFER_WIDTH, FRAMEBUFFER_HEIGHT);
	while (1)
	{
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT) return msg.lParam;
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

        InvalidateRect(window->handle, &windowRect, FALSE);
		Sleep(1);
	}
}

int main()
{
    aligned sse3d_matrix_t matrix;
    aligned sse3d_matrix_t translation;
    aligned sse3d_matrix_t rotation;
    aligned sse3d_matrix_t scale;
    aligned sse3d_matrix_t lookat;

    aligned sse3d_vector_t vector;
    aligned sse3d_vector_t va       = { 32.0f,  1.0f, 1.0f, 0.0f };
    aligned sse3d_vector_t vb       = { 0.0f,   5.5f, 0.5f, 0.0f };
    aligned sse3d_vector_t vc       = { 16.0f,  16.0f, 0.0f, 0.0f };
    aligned sse3d_vector_t target   = { 0.0f, -4.0f, -3.0f, 0.0f };
    aligned sse3d_vector_t camera   = { 0.0f,  0.0f,  1.0f, 0.0f };
    aligned sse3d_vector_t up_vec   = { 0.0f,  1.0f,  0.0f, 0.0f };
    aligned sse3d_vector_t original = { 1.0f,  2.0f,  3.0f, 1.0f };

    int i;
    int nr_matrix_tests = 10000000;
    int nr_vec = 30000000;
    sse3d_vector_t* vector_list;
    clock_t start, end;
    
    int width = 32, height=16;
    unsigned char* scanlines = (unsigned char*) calloc(width * height, 1);

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