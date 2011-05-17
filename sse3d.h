#pragma once

#ifndef SSE3D
#define SSE3D

#include <windows.h>
#include <math.h>
#include <string.h>

#pragma region Types and definitions

#ifndef aligned
#define aligned __declspec(align(16))
#endif

#ifndef aligned_malloc
#define aligned_malloc(size) _aligned_malloc(size, 16)
#endif

#ifndef aligned_free
#define aligned_free(ptr) _aligned_free(ptr)
#endif

#ifndef naked
#define naked __declspec(naked)
#endif

#ifndef FRAMEBUFFER_WIDTH
#define FRAMEBUFFER_WIDTH 800
#endif

#ifndef FRAMEBUFFER_HEIGHT
#define FRAMEBUFFER_HEIGHT 480
#endif

typedef struct
{
    union
    {
        char   i8[16];
        short  i16[8];
        int    i32[4];
        float  f32[4];
        double f64[2];
    };
} sse3d_m128_t;

typedef struct
{
    union
    {
        float _v[4];
        struct
        {
            float x, y, z, w;
        };
    };        
} sse3d_vector_t;

typedef struct
{
    union
    {
        float _v[16];
        struct 
        { 
            float m00, m01, m02, m03,
                  m10, m11, m12, m13,
                  m20, m21, m22, m23,
                  m30, m31, m32, m33;
        };
    };
} sse3d_matrix_t;

typedef void (*pixelshaderproc)(unsigned int *out, sse3d_vector_t *normal);

typedef struct
{
	unsigned short  width, height;
	unsigned float  *z_buffer;
	unsigned int    *p_buffer;
	pixelshaderproc p_shader;
} sse3d_render_ctx_t;

#pragma endregion

#pragma region Forward declarations

void sse3d_multiply_vectors(sse3d_vector_t*, sse3d_matrix_t*, sse3d_vector_t*, int);
void sse3d_multiply_matrix(sse3d_matrix_t*, sse3d_matrix_t*, sse3d_matrix_t*);
void sse3d_transpose_matrix(sse3d_matrix_t*, sse3d_matrix_t*);

void sse3d_normalize_vectors(sse3d_vector_t*, sse3d_vector_t*, int);
float sse3d_dotproduct_vector(sse3d_vector_t*, sse3d_vector_t*);
void sse3d_crossproduct_vector(sse3d_vector_t*, sse3d_vector_t*, sse3d_vector_t*);
void sse3d_prepare_render_vectors(sse3d_vector_t*, int);

void sse3d_identity_matrix(sse3d_matrix_t*);
void sse3d_rotation_x_matrix(sse3d_matrix_t*, float);
void sse3d_rotation_y_matrix(sse3d_matrix_t*, float);
void sse3d_rotation_z_matrix(sse3d_matrix_t*, float);
void sse3d_scale_matrix(sse3d_matrix_t*, float, float, float);
void sse3d_translation_matrix(sse3d_matrix_t*, float, float, float);

#pragma endregion

#pragma region Naked functions

/* ---------------------------------------------- *
 * Naked function to normalize xmm0               *
 * ---------------------------------------------- */
naked _sse3d_normalize_xmm0()
{
    __asm
    {
        movaps  xmm2, xmm0      ;// xmm4 = source vector (x,y,z,0)
        mulps   xmm2, xmm2      ;// xmm4 = (x^2, y^2, z^2, 0)
        haddps  xmm2, xmm2      ;// xmm4 = (x^2 + y^2, z^2 + 0, x^2 + y^2, z^2 + 0)
        haddps  xmm2, xmm2      ;// xmm4 = (x^2 + y^2 + z^2 + 0, ... times 4)
        rsqrtps xmm2, xmm2      ;// xmm4 = (1 / sqrt(x^2 + y^2 + z^2), ... times 4)
        mulps   xmm0, xmm2      ;// xmm0 = for (n in x,y,z,w) n = n/sqrt(x^2 + y^2 + z^2)
        ret
    }
}

/* ---------------------------------------------- *
 * Naked function to calculate the dotproduct     *
 * of xmm0 and xmm1 and store the result in xmm0  *
 * ---------------------------------------------- */
naked _sse3d_dotproduct_xmm0_xmm1()
{
    __asm
    {
        mulps  xmm0, xmm1
        haddps xmm0, xmm0
        haddps xmm0, xmm0
        ret
    }
}


/* ---------------------------------------------- *
 * Naked function to calculate the crossproduct   *
 * of xmm0 and xmm1 and store the result in xmm0  *
 * We calculate the cross product like            *  
 * xmm0.x = xmm0.y * xmm1.z - xmm0.z * xmm1.y,    *
 * xmm0.y = xmm0.z * xmm1.x - xmm0.x * xmm1.z,    *
 * xmm0.z = xmm0.x * xmm1.y - xmm0.y * xmm1.x     *
 * ---------------------------------------------- */
naked _sse3d_crossproduct_xmm0_xmm1()
{
    __asm
    {
        movaps xmm2, xmm0       ;// xmm4 = source vector a
        movaps xmm3, xmm1       ;// xmm5 = source vector b
        
        shufps xmm0, xmm0, 0xC9 ;// 0xC9 = 3,0,2,1 => xmm0 = (a.y, a.z, a.x, a.w)
        shufps xmm2, xmm2, 0xD2 ;// 0xD3 = 3,1,0,2 => xmm4 = (a.z, a.x, a.y, a.w)
        shufps xmm1, xmm1, 0xC9 ;// 0xC9 = 3,0,2,1 => xmm5 = (b.y, b.z, b.x, b.w)
        shufps xmm3, xmm3, 0xD2 ;// 0xD3 = 3,1,0,2 => xmm1 = (b.z, b.x, b.y, b.w)
        
        mulps  xmm0, xmm3       ;// xmm0 = (a.y * b.z, a.z * b.x, a.x * b.y, _w)
        mulps  xmm1, xmm2       ;// xmm1 = (a.z * b.y, a.x * b.z, a.y * b.x, _w)
        subps  xmm0, xmm1       ;// xmm0 = (a.y*b.z - a.z*b.y, 
                                ;//         a.z*b.x - a.x*b.z, 
                                ;//         a.x*b.y - a.y*b.x)
        ret
    }
}

#pragma endregion

#pragma region Vector functions

/* -------------------------------------------------------------------- *
 * Normalize a number of vectors                                        *    
 * Parameters:                                                          *
 * [dst]    A pointer to a memory location holding [nr_vec] vectors     *
 * [src]    A pointer to a memory location holding [nr_vec] vectors     *
 *          which are to be normalized                                  *
 * [nr_vec] The number of vectors in [dst] and [src]                    *
 * -------------------------------------------------------------------- */
void sse3d_normalize_vectors(sse3d_vector_t *dst, sse3d_vector_t *src, int nr_vec)
{
    __asm
    {
        mov esi, src            ;// esi = source vectors
        mov edi, dst            ;// edi = destination vectors
        mov ecx, nr_vec         ;// ecx = number of vectors

        normalize:
        movaps  xmm0, [esi]     ;// xmm0 = current source vector (x,y,z,0)

        call _sse3d_normalize_xmm0;
    
        movaps [edi], xmm0      ;// xmm0 now holds the normalized vector, store in current destination vector
        add esi, 0x10           ;// increment esi to the next vector
        add edi, 0x10           ;// increment edi to the next vector
        loop normalize;         ;// loop while there are vectors
    }
}

/* -------------------------------------------------------------------- *
 * Get the dot product of 2 vectors                                     *
 * -------------------------------------------------------------------- */
float sse3d_dotproduct_vector(sse3d_vector_t *src_a, sse3d_vector_t *src_b)
{
    aligned sse3d_m128_t result;
    
    __asm
    {
        mov esi, src_a          ;// esi = src_a
        mov edi, src_b          ;// edi = src_b
        movaps xmm0, [esi]      ;// xmm0 = src_a
        movaps xmm1, [edi]      ;// xmm1 = src_b
        call _sse3d_dotproduct_xmm0_xmm1;
    
        movaps [result], xmm0;     ;// all f32 values in xmm0 now hold the dotproduct
    }
    return result.f32[0];       // return one of the f32 values
}

/* -------------------------------------------------------------------- *
 * Take the cross product of 2 vectors and store it in dst              *
 * Parameters:                                                          *
 * [dst]    A pointer to a vector that will hold the result             *
 * [src_a]  A pointer to the first vector                               *
 * [src_b]  A pointer to the second vector                              *
 * -------------------------------------------------------------------- */
void sse3d_crossproduct_vector(sse3d_vector_t *dst, sse3d_vector_t *src_a, sse3d_vector_t *src_b)
{
    __asm
    {
        mov esi, src_a          ;// esi = source vector a
        mov edi, src_b          ;// edi = source vector b
        movaps xmm0, [esi]      ;// xmm0 = source vector a (a.x, a.y, a.z, a.w)
        movaps xmm1, [edi]      ;// xmm1 = source vector b (b.x, b.y, b.z, b.w)
        call _sse3d_crossproduct_xmm0_xmm1;
    
        mov edi, dst            ;// edi = destination vector
        movaps [edi], xmm0      ;// xmm0 now holds the cross product, store it in the destination vector
    }
}

/* ---------------------------------------------------------------- *
 * Prepare a number of vectors for rendering by diving by w         *
 * Parameters:                                                      *
 * [vectors]    A pointer to a memory location of [nr_vec] vectors  *
 *              which will be used as input AND output. This is a   *
 *              destructive operation that should be done as the    *
 *              last step before rendering.                         *
 * [nr_vec]     The number of vectors to prepare in [vectors]       *
 * ---------------------------------------------------------------- */
void sse3d_prepare_render_vectors(sse3d_vector_t *vectors, int nr_vec)
{
    __asm
    {
        mov edi, vectors            ;// edi = source and destination vectors
        mov ecx, nr_vec             ;// ecx = number of vectors

        prepare_vector:
        movaps xmm0, [edi]          ;// xmm0 = current vector [x, y, z, w]
        movaps xmm1, xmm0           ;// xmm1 = current vector [x, y, z, w]
        shufps xmm1, xmm1, 0xFF     ;// xmm1 = [w, w, w, w]
        divps  xmm0, xmm1           ;// xmm0 = [x/w, y/w, z/w, 1]
        movaps [edi], xmm0          ;// xmm0 now holds the prepared vector, store it as the current vector
        add edi, 0x10               ;// edi = next vector
        loop prepare_vector         ;// continue while there are vectors
    }
}

#pragma endregion

#pragma region Matrix functions

/* ------------------------------------------------------------------------ *
 * Multiply a number of vectors by a matrix.                                *
 * Parameters:                                                              *
 * [dst]    A pointer to a memory location that can hold [nr_vec] vectors.  *
 * [matrix] A pointer to a matrix used to multiply the [src] vectors.       *
 * [src]    A pointer to a memory location containing the [nr_vec] vectors  *
 *          that need to be multiplied.                                     *
 * [nr_vec] The number of vectors in [src] and [dst].                       *
 * ------------------------------------------------------------------------ */
void sse3d_multiply_vectors(sse3d_vector_t *dst, sse3d_matrix_t *matrix, sse3d_vector_t *src, int nr_vec)
{
    __asm
    {
        mov edi, dst;               ;// edi = destination vectors
        mov esi, matrix             ;// esi = matrix
        movaps xmm0, [esi + 0x00]   ;// xmm0 = matrix row 0
        movaps xmm1, [esi + 0x10]   ;// xmm1 = matrix row 1
        movaps xmm2, [esi + 0x20]   ;// xmm2 = matrix row 2
        movaps xmm3, [esi + 0x30]   ;// xmm3 = matrix row 3

        mov esi, src                ;// esi = source vectors
        mov ecx, nr_vec             ;// ecx = number of vectors in source / destination

        calcvec:
        movaps xmm4, [esi]          ;// xmm4 = current vector
        movaps xmm5, xmm4           ;// xmm5 = current vector
        movaps xmm6, xmm4           ;// xmm6 = current vector
        movaps xmm7, xmm5           ;// xmm7 = current vector
        
        mulps xmm4, xmm0            ;// xmm4 = current vector * matrix row 0
        mulps xmm5, xmm1            ;// xmm5 = current vector * matrix row 1
        mulps xmm6, xmm2            ;// xmm6 = current vector * matrix row 2
        mulps xmm7, xmm3            ;// xmm7 = current vector * matrix row 3

        haddps xmm4, xmm5           ;// horizontal add: xmm4 = (m00 + m01, m02 + m03, m10 + m11, m12 + m13)
        haddps xmm6, xmm7           ;// horizontal add: xmm6 = (m20 + m21, m22 + m23, m30 + m31, m32 + m33)
        haddps xmm4, xmm6           ;// horizontal add: xmm4 = (m00 + m01 + m02 + m03, 
                                    ;//                         m10 + m11 + m12 + m13, 
                                    ;//                         m20 + m21 + m22 + m23, 
                                    ;//                         m30 + m31 + m32 + m33)
        movaps [edi], xmm4          ;// xmm4 now holds the resulting vector, store it in edi
        add edi, 0x10               ;// increment edi to the next vector
        add esi, 0x10               ;// increment esi to the next vector
        loop calcvec                ;// loop if more vectors present
    }
}


/* --------------------------------------------------------------------------------- *
 * Multiply 2 matrices and store the result                                          *
 * Parameters:                                                                       *
 * [dst]    A pointer to a memory location holding the result of the multiplication. *
 * [src_a]  A pointer to the first matrix (row major)                                *
 * [src_b]  A pointer to the second matrix (row major).                              *
 *          This matrix will be transposed into column major and then B*A will be    *
 *          computed giving the row major result of A*B                              *
 * --------------------------------------------------------------------------------- */
void sse3d_multiply_matrix(sse3d_matrix_t *dst, sse3d_matrix_t *src_a, sse3d_matrix_t *src_b)
{
    aligned sse3d_matrix_t transposed;
    sse3d_transpose_matrix(&transposed, src_b);
    sse3d_multiply_vectors((sse3d_vector_t*) dst, &transposed, (sse3d_vector_t*) src_a, 4);
}

/* ----------------------------------------------------------------- *
 * Transpose a matrix from row major to column major or vice versa   *
 * the index mapping is dst_i = (dst_i % 4) * 4 + dst_i / 4          *
 * Parameters:                                                       *
 * [dst] A pointer to a matrix that will hold the transposed matrix. *
 * [src] A pointer to a matrix that is to be transposed.             *
 * ----------------------------------------------------------------- */
void sse3d_transpose_matrix(sse3d_matrix_t *dst, sse3d_matrix_t *src)
{
    __asm
    {
        mov esi, src
        mov edi, dst

        movaps xmm0, [esi + 0x00]   ;// xmm0 = [0 1 2 3]
        movaps xmm1, [esi + 0x10]   ;// xmm1 = [4 5 6 7]
        movaps xmm2, [esi + 0x20]   ;// xmm2 = [8 9 A B]
        movaps xmm3, [esi + 0x30]   ;// xmm3 = [C D E F]

        movaps xmm4, xmm2           ;// xmm4 = [8 9 A B]
        movaps xmm5, xmm2           ;// xmm5 = [8 9 A B]

        punpckldq xmm4, xmm3        ;// xmm4 = [8 C 9 D]
        punpckhdq xmm5, xmm3        ;// xmm5 = [A E B F]
        movaps    xmm2, xmm0        ;// xmm2 = [0 1 2 3]
        punpckldq xmm0, xmm1        ;// xmm0 = [0 4 1 5]
        punpckhdq xmm2, xmm1        ;// xmm2 = [2 6 3 7]
        movaps    xmm1, xmm0        ;// xmm1 = [0 4 1 5]
        movaps    xmm3, xmm2        ;// xmm3 = [2 6 3 7]

        punpcklqdq xmm0, xmm4       ;// xmm0 = [0 4 8 C]
        punpckhqdq xmm1, xmm4       ;// xmm1 = [1 5 9 D]
        punpcklqdq xmm2, xmm5       ;// xmm2 = [2 6 A E]
        punpckhqdq xmm3, xmm5       ;// xmm3 = [3 7 B F]

        movaps [edi + 0x00], xmm0
        movaps [edi + 0x10], xmm1
        movaps [edi + 0x20], xmm2
        movaps [edi + 0x30], xmm3
    }
}

#pragma endregion

#pragma region Matrix transformations

/* ----------------------------------------------- *
 * Load the identity matrix into dst               *
 * ----------------------------------------------- */
void sse3d_identity_matrix(sse3d_matrix_t *dst)
{
    memset(dst, 0, sizeof(float) * 16);
    dst->m00 = 1;
    dst->m11 = 1;
    dst->m22 = 1;
    dst->m33 = 1;
}


/* ---------------------------------------------------- *
 * Load a rotation matrix around the z axis into dst    *
 * Parameters:                                          *
 * [dst]    A pointer to a matrix to load into.         *
 * [angle]  An angle in radians                         *
 * ---------------------------------------------------- */
void sse3d_rotation_z_matrix(sse3d_matrix_t *dst, float angle)
{
    memset(dst, 0, sizeof(float) * 16);
    dst->m00 = (float)cos(angle);
    dst->m01 = (float)-sin(angle);
    dst->m10 = (float)sin(angle);
    dst->m11 = (float)cos(angle);
    dst->m22 = 1.0f;
    dst->m33 = 1.0f;
}


/* ---------------------------------------------------- *
 * Load a rotation matrix around the y axis into dst    *
 * Parameters:                                          *
 * [dst]    A pointer to a matrix to load into.         *
 * [angle]  An angle in radians                         *
 * ---------------------------------------------------- */
void sse3d_rotation_y_matrix(sse3d_matrix_t *dst, float angle)
{
    memset(dst, 0, sizeof(float)*16);
    dst->m00 = (float)cos(angle);
    dst->m02 = (float)-sin(angle);
    dst->m11 = 1.0f;
    dst->m20 = (float)sin(angle);
    dst->m22 = (float)cos(angle);
    dst->m33 = 1.0f;
}


/* ---------------------------------------------------- *
 * Load a rotation matrix around the x axis into dst    *
 * Parameters:                                          *
 * [dst]    A pointer to a matrix to load into.         *
 * [angle]  An angle in radians                         *
 * ---------------------------------------------------- */
void sse3d_rotation_x_matrix(sse3d_matrix_t *dst, float angle)
{
    memset(dst, 0, sizeof(float)*16);
    dst->m00 = 1.0f;
    dst->m11 = (float)cos(angle);
    dst->m12 = (float)-sin(angle);
    dst->m21 = (float)sin(angle);
    dst->m22 = (float)cos(angle);
    dst->m33 = 1.0f;
}


/* ---------------------------------------------------- *
 * Load a scaling matrix into dst                       *
 * Parameters:                                          *
 * [dst]       A pointer to a matrix to load into.      *
 * [scale_(n)] The scale factor for axis (n)            *
 * ---------------------------------------------------- */
void sse3d_scale_matrix(sse3d_matrix_t *dst, float scale_x, float scale_y, float scale_z)
{
    memset(dst, 0, sizeof(float) * 16);
    dst->m00 = scale_x;
    dst->m11 = scale_y;
    dst->m22 = scale_z;
    dst->m33 = 1;
}

/* ---------------------------------------------------- *
 * Load a translation matrix into dst                   *
 * Parameters:                                          *
 * [dst]    A pointer to a matrix to load into.         *
 * [d(n)]   The delta in the (n) axis direction         *
 * ---------------------------------------------------- */
void sse3d_translation_matrix(sse3d_matrix_t *dst, float dx, float dy, float dz)
{
    sse3d_identity_matrix(dst);
    dst->m03 = dx;
    dst->m13 = dy;
    dst->m23 = dz; 
}

#pragma endregion

void sse3d_draw_triangle(sse3d_render_ctx_t* ctx, 
                         sse3d_vector_t *v0, sse3d_vector_t *v1, sse3d_vector_t *v2, 
                         sse3d_vector_t *n0, sse3d_vector_t *n1, sse3d_vector_t *n2)
{
    int i;
    float y;

    aligned sse3d_vector_t  m128_y          = {0, 0, 0, 0};
    aligned sse3d_vector_t  span_from       = {0, 0, 0, 0}, span_to      = {0, 0, 0, 0};
    aligned sse3d_vector_t  normal_from     = {0, 0, 0, 0}, normal_to    = {0, 0, 0, 0};
    aligned sse3d_vector_t  frustrum_min    = {0, 0,-1, 0}, frustrum_max = {FRAMEBUFFER_WIDTH, FRAMEBUFFER_HEIGHT, 1, 0};
    aligned sse3d_vector_t  min, max;

    sse3d_vector_t *va,*vb,*vc,*na,*nb,*nc;

    // Check for degenerate triangles
    if (v0 == v1) return;
    if (v1 == v2) return;
    if (v2 == v0) return;

    // Calculate clipped bounding rect
    __asm
    {
        mov esi, v0
        movaps xmm4, [esi]          ;// xmm4 = v0
        mov esi, v1
        movaps xmm5, [esi]          ;// xmm5 = v1
        mov esi, v2
        movaps xmm6, [esi]          ;// xmm6 = v2

        movaps xmm0, xmm4           ;// xmm0 = v0
        movaps xmm1, xmm4           ;// xmm0 = v0
        minps  xmm0, xmm5           ;// xmm0 = min(v0, v1)
        maxps  xmm1, xmm5           ;// xmm1 = max(v0, v1)
        minps  xmm0, xmm6           ;// xmm0 = min(v0, v1, v2)
        maxps  xmm1, xmm6           ;// xmm1 = max(v0, v1, v2)

        movaps xmm6, [frustrum_min] ;// xmm6 = frustrum min
        movaps xmm7, [frustrum_max] ;// xmm7 = frustrum max
        
        maxps  xmm0, xmm6           ;// xmm0 = max(bounding rect min, frustrum min)
        maxps  xmm1, xmm6           ;// xmm1 = max(bounding rect max, frustrum min)
        minps  xmm0, xmm7           ;// xmm0 = min(bounding rect min, frustrum max)
        minps  xmm1, xmm7           ;// xmm1 = min(bounding rect max, frustrum max)
        
        movaps min, xmm0            ;// min = clipped minima
        movaps max, xmm1            ;// max = clipped maxima   
    }

    // for each scanline from min.y to max.y
    for (y = (float)floor(min.y) + 0.5f; y < max.y; y += 1.0f)
    {
        switch (((v0->y < y) << 2) | ((v1->y < y) << 1) | ((v2->y < y) << 0))
        {
            case 0x0:
            case 0x7: continue;
            case 0x1: va=v1, vb=v2, vc=v0, na=n1, nb=n2, nc=n0; break; 
            case 0x2: va=v0, vb=v1, vc=v2, na=n0, nb=n1, nc=n2; break; 
            case 0x3: va=v1, vb=v0, vc=v2, na=n1, nb=n0, nc=n2; break; 
            case 0x4: va=v2, vb=v0, vc=v1, na=n2, nb=n0, nc=n1; break; 
            case 0x5: va=v2, vb=v1, vc=v0, na=n2, nb=n1, nc=n0; break; 
            case 0x6: va=v0, vb=v2, vc=v1, na=n0, nb=n2, nc=n1; break; 
        }

        m128_y.y = y;
        
        __asm
        {
            mov esi, va
            movaps xmm0, [esi]      ;// xmm0 = 1A
            mov esi, vb
            movaps xmm1, [esi]      ;// xmm1 = 1B
            movaps xmm2, xmm1       ;// xmm2 = 2A
            mov esi, vc
            movaps xmm3, [esi]      ;// xmm3 = 1B
            movaps xmm4, [m128_y]   ;// xmm4 = [0, y, 0, 0]
            movaps xmm6, xmm4       ;// xmm6 = [0, y, 0, 0]
            
            subps  xmm1, xmm0       ;// xmm1 = 1B - 1A = Adelta (Adx, Ady, Adz, 0)
            subps  xmm3, xmm2       ;// xmm3 = 2B - 2A = Bdelta (Bdx, Bdy, Bdz, 0)
            subps  xmm4, xmm0       ;// xmm4 = y - 1A.y
            subps  xmm6, xmm2       ;// xmm6 = y - 2A.y
            divps  xmm4, xmm1       ;// xmm4 = (y - 1A.y)/Ady = Ai (interpolation value)
            divps  xmm6, xmm3       ;// xmm6 = (y - 2A.y)/Bdy = Bi (interpolation value)
            shufps xmm4, xmm4, 0x55 ;// xmm4 = (2, 2, 2, 2) = [Ai, Ai, Ai, Ai]
            shufps xmm6, xmm6, 0x55 ;// xmm6 = (2, 2, 2, 2) = [Bi, Bi, Bi, Bi]
            mulps  xmm1, xmm4       ;// xmm4 = Adelta * Ai
            mulps  xmm3, xmm6       ;// xmm3 = Bdelta * Bi
            addps  xmm0, xmm1       ;// xmm0 = Aintersect (Aix, Aiy, Aiz)
            addps  xmm2, xmm3       ;// xmm2 = Bintersect (Bix, Biy, Biz)
            
            movaps xmm5, [frustrum_min]
            movaps xmm7, [frustrum_max]
            maxps  xmm0, xmm5
            maxps  xmm2, xmm5
            minps  xmm0, xmm7
            minps  xmm2, xmm7

            mov esi, na
            movaps xmm5, [esi]      ;// xmm0 = n1A
            mov esi, nb
            movaps xmm1, [esi]      ;// xmm1 = n1B
            movaps xmm7, xmm1       ;// xmm2 = n2A
            mov esi, nc
            movaps xmm3, [esi]      ;// xmm3 = n1B
            subps  xmm1, xmm5       ;// xmm1 = n1B - n1A = nAdelta (nAdx, nAdy, nAdz, 0)
            subps  xmm3, xmm7       ;// xmm3 = n2B - n2A = nBdelta (nBdx, nBdy, nBdz, 0)
            mulps  xmm1, xmm4       ;// xmm4 = nAdelta * Ai
            mulps  xmm3, xmm6       ;// xmm3 = nBdelta * Bi
            addps  xmm5, xmm1       ;// xmm0 = nAintersect (Aix, Aiy, Aiz)
            addps  xmm7, xmm3       ;// xmm2 = nBintersect (Bix, Biy, Biz)
            
            movaps [span_from], xmm0
            movaps [span_to], xmm2

            movaps xmm0, xmm5
            call _sse3d_normalize_xmm0
            movaps [normal_from], xmm0

            movaps xmm0, xmm7
            call _sse3d_normalize_xmm0
            movaps [normal_to], xmm0
        }

        {
            {
                float *z_buffer_ptr = (float*) (ctx->z_buffer + ((int)y) * ctx->width);
                unsigned char  *n_buffer_ptr = (unsigned char*)  (ctx->p_buffer + ((int)y) * ctx->width);

                float dx = span_to.x - span_from.x;
                float dz = span_to.z - span_from.z;
                aligned sse3d_vector_t currnormal;
                for (i=(int)span_from.x; i<(int)span_to.x; i++)
                {
                    float interpolation = (i - (int)span_from.x) / dx;
                    float depth = interpolation * dz + span_from.z;
                    
                    if (depth>1) continue;
                    if (depth<-1) continue;
                    
                    depth+= 1;
                    depth/=2;
                    
                    if (depth > z_buffer_ptr[i]) 
                    {
						aligned sse3d_vector_t light = {0, -0.1f, 1, 0};
						static aligned sse3d_matrix_t rotation;
                        float r, g, b, v;
						currnormal.x = (normal_from.x + (normal_to.x - normal_from.x) * interpolation);
                        currnormal.y = (normal_from.y + (normal_to.y - normal_from.y) * interpolation);
                        currnormal.z = (normal_from.z + (normal_to.z - normal_from.z) * interpolation);
						currnormal.w = 0;
                        sse3d_normalize_vectors(&currnormal,&currnormal, 1);
						sse3d_normalize_vectors(&light, &light,1);
						v = sse3d_dotproduct_vector(&currnormal, &light); //sqrt(currnormal.z + currnormal.y);
						
						if (v < 0) v = 0;
						v = (pow(v, 6));
						v = (sin(v*M_PI) + 5) / 6;
                        r = ((1 + currnormal.y) / 2) * v;
                        g = ((1 + currnormal.x) / 2) * v;
                        b = ((1 + currnormal.z) / 2) * v;

                        z_buffer_ptr[i] = depth;
                        
                        n_buffer_ptr[i*4+0] = ((unsigned char)(r * 0xFF * depth));
                        n_buffer_ptr[i*4+1] = ((unsigned char)(g * 0xFF * depth));
                        n_buffer_ptr[i*4+2] = ((unsigned char)(b * 0xFF * depth));
                        
                    }
                }
            }
        }
    }
}

#endif