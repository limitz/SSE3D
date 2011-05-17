#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "sse3d.h"
#include "sse3d_window.h"
#include "sse3d_dragon.h"

static aligned sse3d_vector_t v_transformed[dragon_nr_vertices];
static aligned sse3d_vector_t n_transformed[dragon_nr_normals];
static aligned sse3d_vector_t light_direction = {0, 0.0f, 1.0f, 0};

void pixel_shader_depth(unsigned char *color, sse3d_vector_t *normal, float depth)
{
    color[0] = color[1] = color[2] =(unsigned char) (depth * 0xFF);
}

void pixel_shader_normals(unsigned char *color, sse3d_vector_t *normal, float depth)
{
    color[0] = (unsigned char) ((normal->x + 1) * 0x7F);
    color[1] = (unsigned char) ((normal->y + 1) * 0x7F);
    color[2] = (unsigned char) ((normal->z + 1) * 0x7F);
}

void pixel_shader_dark(unsigned char *color, sse3d_vector_t *normal, float depth)
{
    float v = sse3d_dotproduct_vector(&light_direction, normal);
    if (v < 0) v = 0;
    color[0] = color[1] = color[2] =(unsigned char) (pow(v,4.0f) * depth * 0xFF);;
}

void pixel_shader_chrome(unsigned char *color, sse3d_vector_t *normal, float depth)
{
    float v = sse3d_dotproduct_vector(&light_direction, normal);
    if (v < 0) v = (-v);
    color[0] = color[1] = color[2] =(unsigned char) (pow(v,4.0f) * depth * 0x300);
}

void pixel_shader_colors(unsigned char *color, sse3d_vector_t *normal, float depth)
{
    float v = sse3d_dotproduct_vector(&light_direction, normal);
    if (v < 0) v = 0;
    
    color[0] = (unsigned char) ((normal->z + 1) * pow(v,3) * depth * 0x7F);
    color[1] = (unsigned char) ((normal->x + 1) * pow(v,3) * depth * 0x7F);
    color[2] = (unsigned char) ((normal->y + 1) * pow(v,3) * depth * 0x7F);
}

void pixel_shader_silk(unsigned char *color, sse3d_vector_t *normal, float depth)
{
    float v = sse3d_dotproduct_vector(&light_direction, normal);
    float r;
    

    r = (sin(normal->x*7) + cos(normal->y*9) + sin(depth*5)) / 6 + 0.5;
    
    color[0] = (unsigned char) (r * pow(v,2) * depth * 0x7F);
    color[1] = (unsigned char) (r * pow(v,2) * depth * 0x7F);
    color[2] = (unsigned char) (r * pow(v,2) * depth * 0x7F);
}

void render(unsigned float *z_buffer, unsigned int *n_buffer, int width, int height)
{
    
    int i;
    static float angle = 0;
    static aligned sse3d_matrix_t light_rotation;
    static aligned sse3d_matrix_t model, model_scale, model_rotation, model_rotation_y, model_rotation_x,model_translation;
    static aligned sse3d_matrix_t projection, projection_scale, projection_translation, lookat;
    static aligned sse3d_matrix_t identity, transform;

    sse3d_rotation_y_matrix(&light_rotation, -0.02f);
    //sse3d_multiply_vectors(&light_direction, &light_rotation, &light_direction, 1);
    sse3d_normalize_vectors(&light_direction, &light_direction, 1);
    
    sse3d_identity_matrix(&identity);
    sse3d_scale_matrix(&model_scale, 1, 1, 1);
    sse3d_rotation_x_matrix(&model_rotation, -M_PI/2.0f);
    sse3d_rotation_x_matrix(&model_rotation_x, 0);
    sse3d_rotation_y_matrix(&model_rotation_y, angle += 0.005f);
    sse3d_translation_matrix(&model_translation, 0.0f, -0.1f, 0.15f);
    
    sse3d_multiply_matrix(&model, &model_rotation_y, &model_rotation);
    sse3d_multiply_matrix(&model, &model_rotation_x, &model);
    sse3d_multiply_matrix(&model, &model_scale, &model);
    sse3d_multiply_matrix(&model, &model_translation, &model);
    
    sse3d_scale_matrix(&projection_scale, 3.7f * height / 2.0f, 3.7f * height / 2.0f, 1.0f);
    sse3d_translation_matrix(&projection_translation, width / 2.0f, height / 2.0f, 0.0f);
    sse3d_multiply_matrix(&projection, &projection_translation, &projection_scale);

    sse3d_multiply_matrix(&transform, &projection, &model);
    sse3d_multiply_vectors(v_transformed, &transform, dragon_vertices, dragon_nr_vertices);
    sse3d_multiply_vectors(n_transformed, &model, dragon_normals, dragon_nr_normals);

    sse3d_prepare_render_vectors(v_transformed, dragon_nr_vertices);
    for (i=0; i<dragon_nr_indices; i+=6)
    {
        sse3d_render_ctx_t ctx;
        ctx.width = width;
        ctx.height = height;
        ctx.p_buffer = n_buffer;
        ctx.z_buffer = z_buffer;
        ctx.p_shader = 0;
        ctx.p_shader = pixel_shader_silk;

        sse3d_draw_triangle(&ctx, 
            v_transformed + dragon_indices[i+0], v_transformed + dragon_indices[i+1], v_transformed + dragon_indices[i+2],
            n_transformed + dragon_indices[i+3], n_transformed + dragon_indices[i+4], n_transformed + dragon_indices[i+5]);
    }
}

int swapat = 0;
int WINAPI WinMain(HINSTANCE instance, HINSTANCE prevInstance, LPSTR cmdLine, int showCmd)
{
    int i;
    MSG msg;
    RECT windowRect;
    char buffer[1024];

    sse3d_window_t* window = sse3d_create_window(instance, render);
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