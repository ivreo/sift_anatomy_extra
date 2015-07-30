#ifndef _LIB_UTIL_H_
#define _LIB_UTIL_H_ 

#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT2
    #define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_LN2
    #define M_LN2 0.693147180559945309417
#endif

//#define QUAD 1
#define DOUBLE 1


// FLOATING POINT PRECISION
#ifdef QUAD
    #include <quadmath.h>
    typedef __float128 _myfloat;
#else
    #ifdef DOUBLE
        typedef double _myfloat;
    #else
        typedef float _myfloat;
    #endif
#endif





void* xmalloc(size_t size);

// Reallocate memory of abort on failure
void* xrealloc(void* p, size_t size);

// Free memory allocated by xmalloc or xrealloc.
void xfree(void* p);

// Write a formatted message to standard error and abort the program, for
// example, fatal_error("Failed to open file %s for writing", filename);
void fatal_error(const char* format, ...);

// Write formatted message to standard error for debugging.
void debug(const char* format, ...);

// Linearly rescale input such that its values are in the range [0, 250].
void linear_conversion(const _myfloat *in, _myfloat *out, int l);

// Find the maximum value of an array.
_myfloat array_max(const _myfloat* array, int length);

// Find the minimum value of an array.
_myfloat array_min(const _myfloat* array, int length);

// Find the maximum value of an array.
_myfloat find_array_max(const _myfloat* array, int length, int *position);

// Find the minimum value of an array.
_myfloat find_array_min(const _myfloat* array, int length, int *position);

// Find the two minimal values in an array.
void find_array_two_min(const _myfloat* array, int length, _myfloat* minA, _myfloat* minB, int* iA, int* iB);

// L2 norm of an array.
_myfloat array_l2norm(const _myfloat* array, int length);

// Compute the SQUARE of the euclidean distance
_myfloat euclidean_distance_square(const _myfloat* a1, const _myfloat* a2, int length);

// Compute the euclidean distance
_myfloat euclidean_distance(const _myfloat* a1, const _myfloat* a2, int length);

// Compute the x modulus y
_myfloat modulus(_myfloat x, _myfloat y);

// Multiply the rotation matric R_alpha with [x,y]^T
void apply_rotation(_myfloat x, _myfloat y, _myfloat *rx, _myfloat *ry, _myfloat alpha);

#endif // _LIB_UTIL_H_
