#pragma once

#ifndef SSE3D_WINDOW
#define SSE3D_WINDOW

#ifndef FRAMEBUFFER_WIDTH
#define FRAMEBUFFER_WIDTH 480
#endif

#ifndef FRAMEBUFFER_HEIGHT
#define FRAMEBUFFER_HEIGHT 320
#endif

#include <windows.h>
#include "sse3d.h"

typedef void (*renderproc)(unsigned char *buffer, int width, int height);

typedef struct
{
    HWND handle;
    HDC backbuffer_dc;
    HBITMAP backbuffer;
    unsigned char* buffer;
    renderproc render;
} sse3d_window_t;



LRESULT CALLBACK sse3d_windowproc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    sse3d_window_t* instance = (sse3d_window_t*) GetWindowLong(hwnd, GWL_USERDATA);

	switch (msg)
	{
		case WM_CREATE:
        {
            CREATESTRUCT* create = (CREATESTRUCT*) lParam;
			instance = (sse3d_window_t*) create->lpCreateParams;
            instance->handle = hwnd;
            
			SetWindowLong(hwnd, GWL_USERDATA, (LONG)instance);
			SetFocus(hwnd);
            return 0;
        }
        case WM_CLOSE:
		{
			DestroyWindow(hwnd);
			PostQuitMessage(0);
			return 0;
		}
		
        case WM_DESTROY: 
			PostQuitMessage(0); 
			return 0;
	
        case WM_ERASEBKGND:
			return 1;

		case WM_PAINT:
			if (GetUpdateRect(hwnd, NULL, FALSE))
			{
				PAINTSTRUCT paint;
				HDC hdc = BeginPaint(hwnd, &paint);

				RECT invalidRect = paint.rcPaint;
				SetViewportOrgEx(hdc, 0, 0, NULL);
                memset(instance->buffer, 0x00, FRAMEBUFFER_WIDTH * FRAMEBUFFER_HEIGHT);
                if (instance->render) instance->render(instance->buffer, FRAMEBUFFER_WIDTH, FRAMEBUFFER_HEIGHT);

                BitBlt(hdc, 0, 0, FRAMEBUFFER_WIDTH, FRAMEBUFFER_HEIGHT, instance->backbuffer_dc, 0, 0, SRCCOPY);
				EndPaint(hwnd, &paint);
			}
			return 0;

		default:
			return DefWindowProc(hwnd, msg, wParam, lParam);
	}
}


sse3d_window_t* sse3d_create_window(HINSTANCE instance, renderproc render)
{
    int i;
            
    HDC hdc;
    BITMAPINFO* bminfo;
    RGBQUAD* palette;
	RECT windowRect;

    WNDCLASS windowClass;
    sse3d_window_t* result = (sse3d_window_t*) calloc(1, sizeof(sse3d_window_t));

    memset(&windowClass, 0, sizeof(WNDCLASS));
	windowClass.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	windowClass.lpfnWndProc = sse3d_windowproc;
	windowClass.hInstance = instance;
	windowClass.lpszClassName = "SSE3D";
    windowClass.hCursor = LoadCursor(NULL, IDC_ARROW);

    if (!RegisterClass(&windowClass)) return 0;

    SetRect(&windowRect, 0, 0, FRAMEBUFFER_WIDTH, FRAMEBUFFER_HEIGHT);
    AdjustWindowRect(&windowRect,WS_POPUP | WS_CAPTION | WS_SYSMENU, 0);

    CreateWindow(windowClass.lpszClassName, windowClass.lpszClassName, 
                 WS_POPUP | WS_CAPTION | WS_SYSMENU, 
                 CW_USEDEFAULT, CW_USEDEFAULT, windowRect.right - windowRect.left, windowRect.bottom - windowRect.top,
                 NULL, NULL, windowClass.hInstance, result);

    if (!result->handle)
    {
        free(result);
        return 0;
    }

    hdc = GetDC(result->handle);
    result->backbuffer_dc = CreateCompatibleDC(hdc);
    result->render = render;
    bminfo = (BITMAPINFO*) calloc(1, sizeof(BITMAPINFO) + sizeof(RGBQUAD) * 256);
    bminfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bminfo->bmiHeader.biBitCount = 8;
    bminfo->bmiHeader.biWidth =  FRAMEBUFFER_WIDTH;
    bminfo->bmiHeader.biHeight = FRAMEBUFFER_HEIGHT;
    bminfo->bmiHeader.biPlanes = 1;
    bminfo->bmiHeader.biCompression = BI_RGB;

    palette = (RGBQUAD*)(((char*)bminfo) + bminfo->bmiHeader.biSize);

    for (i=0; i<256; i++)
    {
        palette[i].rgbRed = i;
        palette[i].rgbGreen = i;
        palette[i].rgbBlue = i;
    }

    result->backbuffer = CreateDIBSection(result->backbuffer_dc, bminfo, DIB_RGB_COLORS, &result->buffer, NULL, NULL);
    SelectObject(result->backbuffer_dc, result->backbuffer);
    free(bminfo);
    ReleaseDC(result->handle, hdc);

    ShowWindow(result->handle, SW_NORMAL);
    UpdateWindow(result->handle);
    return result;
}




#endif
