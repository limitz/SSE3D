/* ------------------------------------------ *
 * SSE3D - 3D renderer using SSE instructions *
 * ------------------------------------------ *
 * Author: Eddy Pronk                         *
 * Date:   18 May 2011                        *
 * ------------------------------------------ */

#define WIDTH 800
#define HEIGHT 600

#pragma region Includes

#pragma comment(lib, "msimg32")

#include <windows.h>
#include <windowsx.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "sse3d.h"
#include "sse3d_dragon.h"

#pragma endregion

#pragma region Static fields

static aligned sse3d_render_params_t s_params;
static aligned sse3d_model_t s_model;

static HDC		s_hdc;
static HBITMAP	s_hbm;

static float    s_delta_angle = 0.002f;
static float    s_speed = 4.0f;
static float	s_alpha = 0.5f;
static float	s_correct = 2.0f;

#pragma endregion

#pragma region Pixel shaders

/* -------------------------------------------- *
 * PS: Draw the z-buffer values                 *
 * -------------------------------------------- */
void pixel_shader_depth(sse3d_render_params_t* params, unsigned char *color, sse3d_vector_t *normal, float depth)
{
    color[0] = color[1] = color[2] =(unsigned char) (depth * 0xFF);
}

/* -------------------------------------------- *
 * PS: Draw the raw normal values               *
 * -------------------------------------------- */
void pixel_shader_normals(sse3d_render_params_t* params, unsigned char *color, sse3d_vector_t *normal, float depth)
{
    color[0] = (unsigned char) ((normal->z + 1.0f) * 0x7F);
    color[1] = (unsigned char) ((normal->y + 1.0f) * 0x7F);
    color[2] = (unsigned char) ((normal->x + 1.0f) * 0x7F);
}

/* ------------------------------------------------ *
 * PS: Draw the model with light intensity and fog  *
 * ------------------------------------------------ */
void pixel_shader_dark(sse3d_render_params_t* params, unsigned char *color, sse3d_vector_t *normal, float depth)
{
    float v = sse3d_dotproduct_vector(&params->light_direction, normal);
    if (v < 0) v = 0;
    color[0] = color[1] = color[2] =(unsigned char) (pow(v,4.0f) * depth * 0xFF);;
}

/* -------------------------------------------- *
 * PS: Draw the model using normal data         *
 * -------------------------------------------- */
void pixel_shader_colors(sse3d_render_params_t* params, unsigned char *color, sse3d_vector_t *normal, float depth)
{
    float v = sse3d_dotproduct_vector(&params->light_direction, normal);
    if (v < 0) v = 0;
    
    color[0] = (unsigned char) ((normal->z + 1.0f) * pow(v,3.0f) * depth * 0x7F);
    color[1] = (unsigned char) ((normal->x + 1.0f) * pow(v,3.0f) * depth * 0x7F);
    color[2] = (unsigned char) ((normal->y + 1.0f) * pow(v,3.0f) * depth * 0x7F);
}

/* ----------------------------------------------- *
 * PS: Draw the model using normal data and gloss  *
 * ----------------------------------------------- */
void pixel_shader_silk(sse3d_render_params_t* params, unsigned char *color, sse3d_vector_t *normal, float depth)
{
    float v = sse3d_dotproduct_vector(&params->light_direction, normal);
    float r;

    r = (float)((sin(normal->x*7.0f) + cos(normal->y*9.0f) + sin(params->light_direction.x*5.0f) + sin(depth*12.0f)) / 8.0f + 0.5f);
    
    color[0] = (unsigned char) (r * pow(v,4.0f) * depth * (1.0f+normal->x)/4.0f * 0xFF);
    color[1] = (unsigned char) (r * pow(v,2.0f) * depth * (1.0f+normal->y)/2.0f * 0xFF);
    color[2] = (unsigned char) (r * pow(v,2.0f) * depth * 0xFF);
}

#pragma endregion


/* ------------------------------------------ *
 * Render the dragon model                    *
 * Parameters:                                *
 * [dest] -> Destination DC to draw to        *
 * ------------------------------------------ */
void render(HDC dest)
{
	BLENDFUNCTION bf;
    static float angle = 0;
    static aligned sse3d_matrix_t light_rotation;
    static aligned sse3d_matrix_t model, model_scale, model_rotation, model_rotation_y, model_rotation_x,model_translation;
    static aligned sse3d_matrix_t projection, projection_scale, projection_translation, lookat;
    static aligned sse3d_matrix_t identity, transform;

	sse3d_rotation_x_matrix(&model_rotation_x, -M_PI/2.0f);
    sse3d_rotation_y_matrix(&model_rotation_y, angle += s_delta_angle * s_correct * s_speed);
	sse3d_translation_matrix(&model_translation, 0.0f, -0.1f, 0.15f);
    
    sse3d_multiply_matrix(&s_params.model_matrix, &model_rotation_y, &model_rotation_x);
    sse3d_multiply_matrix(&s_params.model_matrix, &model_translation, &s_params.model_matrix);

	sse3d_clear_buffers(&s_params);
    sse3d_render_model(&s_params, &s_model);

	bf.BlendOp = AC_SRC_OVER;
	bf.BlendFlags = 0;
	bf.SourceConstantAlpha = (unsigned char)(s_alpha * 0xFF);
	bf.AlphaFormat = 0;

    AlphaBlend(dest, 0, 0, s_params.width, s_params.height, s_hdc, 0, 0, s_params.width, s_params.height, bf);
}


#pragma region Window message process

LRESULT CALLBACK window_msg_proc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    static int isMouseDown = 0;

	switch (msg)
	{
        case WM_CLOSE:
			DestroyWindow(hwnd);
        case WM_DESTROY: 
			PostQuitMessage(0); 
			return 0;
	
		case WM_KEYDOWN:
			switch (wParam)
	        {
		        case '0': s_params.p_shader = 0;					s_correct = 0.7f; break;
		        case '1': s_params.p_shader = pixel_shader_depth;	s_correct = 0.7f; break;
		        case '2': s_params.p_shader = pixel_shader_normals; s_correct = 1.0f; break;
		        case '3': s_params.p_shader = pixel_shader_dark;	s_correct = 1.2f; break;
		        case '4': s_params.p_shader = pixel_shader_colors;  s_correct = 1.8f; break;
		        case '5': s_params.p_shader = pixel_shader_silk;	s_correct = 2.0f; break;

		        case ' ': s_delta_angle = s_delta_angle > 0 ? 0 : 0.002f; break;
		        case 'R': s_speed = -s_speed; break;
		        case 'Z': s_alpha -= 0.05f; if (s_alpha < 0.15f) s_alpha = 0.15f; break;
		        case 'X': s_alpha += 0.05f; if (s_alpha > 1.0f) s_alpha = 1.0f; break;

		        case VK_OEM_MINUS:
		        case VK_SUBTRACT: s_speed /= 1.2f; break;

		        case VK_OEM_PLUS:
		        case VK_ADD: s_speed *= 1.2f; break;
            }
            return 0;
			
		case WM_LBUTTONDOWN:
				isMouseDown = 1;
				SetCapture(hwnd);
				return 0;

		case WM_LBUTTONUP:
				isMouseDown = 0;
				ReleaseCapture();
				return 0;

		case WM_MOUSEMOVE:
				if (isMouseDown) 
                {
                    int x = GET_X_LPARAM(lParam);
                    int y = GET_Y_LPARAM(lParam);
                    s_params.light_position.x = 1 - (x / (float)s_params.width)*2;
	                s_params.light_position.y = (y / (float)s_params.height)*2 - 1;
	                s_params.light_position.z = -1.0f;
	                __asm
	                {
		                mov esi, s_params
		                movaps xmm1, [s_params + 0x00];
		                movaps xmm0, [s_params + 0x10];
		                subps  xmm0, xmm1;
		                call _sse3d_normalize_xmm0
		                movaps [s_params+0x20], xmm0;
	                }
                }
				return 0;

        case WM_ERASEBKGND:
			return 1;

		case WM_PAINT:
			if (GetUpdateRect(hwnd, NULL, FALSE))
			{
				PAINTSTRUCT paint;
				HDC hdc = BeginPaint(hwnd, &paint);
				render(hdc);
				EndPaint(hwnd, &paint);
			}
			return 0;

		default:
			return DefWindowProc(hwnd, msg, wParam, lParam);
	}
}

#pragma endregion

#pragma region Initialization

/* -------------------------------------------- *
 * Initialize the render parameters				*
 * Parameters:									*
 * [window] -> The window handle of the screen	*
 * -------------------------------------------- */
int initialize_renderer(HWND window, int width, int height)
{
	HDC windowDC = GetDC(window);
	
	aligned sse3d_matrix_t projection_scale;
	aligned sse3d_matrix_t projection_translation;

	BITMAPINFO bminfo = {0};
	bminfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	bminfo.bmiHeader.biWidth = width;
	bminfo.bmiHeader.biHeight = height;
	bminfo.bmiHeader.biBitCount = 32;
    bminfo.bmiHeader.biPlanes = 1;
    bminfo.bmiHeader.biCompression = BI_RGB;

	s_hdc = CreateCompatibleDC(windowDC);
    if (!s_hdc) return -1;

	s_hbm = CreateDIBSection(s_hdc, &bminfo, DIB_RGB_COLORS, &s_params.p_buffer, NULL, 0);
    if (!s_hbm) return -2;

	SelectObject(s_hdc, s_hbm);

    s_params.width = width;
	s_params.height = height;
	s_params.z_buffer = (float*) calloc(width * height, sizeof(float));
    if (!s_params.z_buffer) return -3;
	
	ReleaseDC(window, windowDC);

	s_params.light_position.z = -1;
	s_params.light_direction.z = 1;

	sse3d_identity_matrix(&s_params.model_matrix);
	sse3d_identity_matrix(&s_params.projection_matrix);

	sse3d_scale_matrix(&projection_scale, 2.5f * height/2.0f, 2.5f * height/2.0f, 1.0f);
	sse3d_translation_matrix(&projection_translation, width/2.0f, height/2.0f, 0);
	sse3d_multiply_matrix(&s_params.projection_matrix, &projection_translation, &projection_scale);

	s_params.p_shader = pixel_shader_silk;

	s_model.nr_vertices = dragon_nr_vertices;
	s_model.nr_normals = dragon_nr_normals;
	s_model.nr_indices = dragon_nr_indices;

	s_model.vertices = dragon_vertices;
	s_model.normals = dragon_normals;
	s_model.indices = dragon_indices;

    return 0;
}

/* ----------------------------------------------------- *
 * Initialize the window                                 *
 * Parameters:									         *
 * [instance] -> HINSTANCE of the application            *
 * [width]    -> Width of client size of the window      * 
 * [height]   -> Height of client size of the window     * 
 * ----------------------------------------------------- */
HWND initialize_window(HINSTANCE instance, int width, int height)
{
    RECT windowRect;

    WNDCLASS windowClass = {0};
	windowClass.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	windowClass.lpfnWndProc = window_msg_proc;
	windowClass.hInstance = instance;
	windowClass.lpszClassName = "SSE3D";
    windowClass.hCursor = LoadCursor(NULL, IDC_ARROW);

    if (!RegisterClass(&windowClass)) return 0;

    SetRect(&windowRect, 0, 0, width, height);
    AdjustWindowRect(&windowRect, WS_POPUP | WS_CAPTION | WS_SYSMENU, 0);

    return CreateWindow(windowClass.lpszClassName,  "SSE3D  --  [0-5] = Pixel Shader,  [+/-] = Speed,  [Z/X] = Alpha,  [Space] = Stop,  [R] = Reverse,  [Mouse] = Light Position", 
                 WS_POPUP | WS_CAPTION | WS_SYSMENU, 
                 CW_USEDEFAULT, CW_USEDEFAULT, windowRect.right - windowRect.left, windowRect.bottom - windowRect.top,
                 NULL, NULL, windowClass.hInstance, NULL);
}

#pragma endregion


/* -------------------------------------------- *
 * Main entry point of the application          *
 * -------------------------------------------- */
int WINAPI WinMain(HINSTANCE instance, HINSTANCE prevInstance, LPSTR cmdLine, int showCmd)
{
    int error;

    MSG msg;
    RECT invalidRect;
    HWND hwnd = initialize_window(instance, WIDTH, HEIGHT);
    if (!hwnd) return 1;

	error = initialize_renderer(hwnd, WIDTH, HEIGHT);
    if (error) return error;
	
    GetClientRect(hwnd, &invalidRect);
    ShowWindow(hwnd, SW_NORMAL);
    UpdateWindow(hwnd);
    
	while (1)
	{
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT) return msg.lParam;
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
        InvalidateRect(hwnd, &invalidRect, FALSE);
		Sleep(1);
	}
}