#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "sse3d.h"
#include "sse3d_window.h"

#include "sse3d_dragon.h"

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
aligned int indices[30000];
int nr_vertices, nr_normals, nr_indices;

aligned sse3d_vector_t v_transformed[10000];
aligned sse3d_vector_t n_transformed[10000];

void render(unsigned float *z_buffer, unsigned int *n_buffer, int width, int height)
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
    sse3d_rotation_y_matrix(&model_rotation_y, angle += 0.0016f);
    sse3d_translation_matrix(&model_translation, 0.0f, -0.1f, 0.15f);
    
    sse3d_multiply_matrix(&model, &model_rotation_y, &model_rotation);
    sse3d_multiply_matrix(&model, &model_rotation_x, &model);
    sse3d_multiply_matrix(&model, &model_scale, &model);
    sse3d_multiply_matrix(&model, &model_translation, &model);
    
    sse3d_scale_matrix(&projection_scale, 3.7f * height / 2.0f, 3.7f * height / 2.0f, 1.0f);
    sse3d_translation_matrix(&projection_translation, width / 2.0f, height / 2.0f, 0.0f);
    sse3d_multiply_matrix(&projection, &projection_translation, &projection_scale);

    sse3d_multiply_matrix(&transform, &projection, &model);
    sse3d_multiply_vectors(v_transformed, &transform, vertices, nr_vertices);
    sse3d_multiply_vectors(n_transformed, &model, normals, nr_normals);

    sse3d_prepare_render_vectors(v_transformed, nr_vertices);
    for (i=0; i<nr_indices; i+=6)
    {
        sse3d_render_ctx_t ctx;
        ctx.width = width;
        ctx.height = height;
        ctx.p_buffer = n_buffer;
        ctx.z_buffer = z_buffer;

        sse3d_draw_triangle(&ctx, 
            v_transformed + indices[i+0] - 1, v_transformed + indices[i+1] - 1, v_transformed + indices[i+2] - 1,
            n_transformed + indices[i+3] - 1, n_transformed + indices[i+4] - 1, n_transformed + indices[i+5] - 1);
    }
}

int WINAPI WinMain(HINSTANCE instance, HINSTANCE prevInstance, LPSTR cmdLine, int showCmd)
{
    int i;
    MSG msg;
    RECT windowRect;
    char buffer[1024];

    sse3d_window_t* window;

    
    FILE *dragon = fopen("dragon.txt", "r");
    FILE* output = fopen("sse3d_dragon.h", "w");
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
    fprintf(output, "#pragma once\n\n#ifndef SSE3D_DRAGON\n#define SSE3D_DRAGON\n\n");
    fprintf(output, "const int dragon_nr_vertices = %d;\n",  nr_vertices);
    fprintf(output, "const int dragon_nr_normals = %d;\n",  nr_normals);
    fprintf(output, "const int dragon_nr_indices = %d;\n\n",  nr_indices);

   
    fprintf(output, "aligned sse3d_vector_t dragon_vertices[] = {\n    ");
    for (i=0; i<nr_vertices; i++)
    {
        fprintf(output, "{%s%0.12ff,%s%0.12ff,%s%0.12ff, 1.0f }, ", vertices[i].x > 0 ? " " : "", vertices[i].x *  .02f, vertices[i].y > 0 ? " " : "", vertices[i].y *  .02f, vertices[i].z > 0 ? " " : "", vertices[i].z *  .02f);
        if (i%8 == 7) fprintf(output, "\n    ");
    }
    fprintf(output, "};\n\n");
    
    fprintf(output, "aligned sse3d_vector_t dragon_normals[] = {\n    ");
    for (i=0; i<nr_normals; i++)
    {
        fprintf(output, "{%s%0.12ff,%s%0.12ff,%s%0.12ff, 0.0f }, ", normals[i].x > 0 ? " " : "", normals[i].x, normals[i].y > 0 ? " " : "", normals[i].y, normals[i].z > 0 ? " " : "", normals[i].z);
        if (i%8 == 7) fprintf(output, "\n    ");
    }
    fprintf(output, "};\n\n");

    fprintf(output, "aligned sse3d_vector_t dragon_indices[] = {\n    ");
    for (i=0; i<nr_indices; i++)
    {
        fprintf(output, "%5d, ", indices[i]);
        if (i%24 == 23) fprintf(output, "\n    ");
    }
    fprintf(output, "};\n\n#endif\n\n");


    fclose(output);
    
    
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
		//Sleep(10);
	}
}