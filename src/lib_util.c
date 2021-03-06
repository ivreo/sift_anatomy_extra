#include <float.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib_util.h"

// Allocate memory or abord on failure
void* xmalloc(size_t size)
{
    if (size == 0)
        fprintf(stderr,"xmalloc: zero size");
    void *p = malloc(size);
    if (!p)
    {
        double sm = size / (0x100000 * 1.0);
        fprintf(stderr,"xmalloc: out of memory when requesting "
                "%zu bytes (%gMB)",//:\"%s\"",
                size, sm);//, strerror(errno));
    }
    return p;
}

// Reallocate memory of abort on failure
void* xrealloc(void* p, size_t size)
{
    void *r = realloc(p, size);
    if (!r) fprintf(stderr,"realloc failed");
    return r;
}

// Free memory allocated by xmalloc or xrealloc.
void xfree(void* p)
{
    if (!p)
        fprintf(stderr,"trying to free a null pointer");
    free(p);
}

// Write a formatted message to standard error and abort the program, for
// example, fatal_error("Failed to open file %s for writing", filename);
void fatal_error(const char* format, ...)
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "Fatal error: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
    exit(1);
}

// Write formatted message to standard error for debugging.
void debug(const char* format, ...)
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "Debug: ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

// Find the maximum value of an array and store its index in the array.
_myfloat find_array_max(const _myfloat* array, int length, int *position)
{
    double max = -FLT_MAX;
    for(int i = 0; i < length; i++){
        if(array[i] > max){
            *position = i;
            max = array[i];
        }
    }
    return max;
}

// Find the minimum value of an array and store its index in the array.
_myfloat find_array_min(const _myfloat* array, int length, int *position)
{
    double min = FLT_MAX;
    for(int i = 0; i < length; i++){
        if(array[i] < min){
            min = array[i];
            *position = i;
        }
    }
    return min;
}

// Find the maximum value of an array.
_myfloat array_max(const _myfloat* array, int length)
{
    int i;
    _myfloat max = find_array_max(array, length, &i);
    (void)i;
    return max;
}

// Find the minimum value of an array.
_myfloat array_min(const _myfloat* array, int length)
{
    int i;
    _myfloat min = find_array_min(array, length, &i);
    (void)i;
    return min;
}

// Find the two minimal values of an array.
void find_array_two_min(const _myfloat* array, int length, _myfloat* minA, _myfloat* minB, int* iA, int* iB)
{
    if(array[0] < array[1]){
        *iA = 0;
        *iB = 1;
        *minA = array[0];
        *minB = array[1];
    }
    else{
        *iA = 1;
        *iB = 0;
        *minA = array[1];
        *minB = array[0];
    }
    for(int i = 2; i < length; i++){
        if(array[i] < *minA){
            *minB = *minA;
            *minA = array[i];
            *iB = *iA;
            *iA = i;
        }
        else if( array[i] < *minB){
            *minB = array[i];
            *iB = i;
        }
    }
}

// Compute the SQUARE Euclidean distance.
_myfloat euclidean_distance_square(const _myfloat* a1, const _myfloat* a2, int length)
{
    _myfloat d = 0.0;
    for(int i = 0; i < length; i++){
        _myfloat t = (a1[i] - a2[i]);
        d += t*t;
    }
    return d;
}

// Compute the Euclidean distance.
_myfloat euclidean_distance(const _myfloat* a1, const _myfloat* a2, int length)
{
    _myfloat d = euclidean_distance_square(a1,a2,length);
    d = sqrt(d);
    return d;
}

// L2 norm of a vector
_myfloat array_l2norm(const _myfloat* array, int length)
{
    _myfloat l2norm = 0;
    for(int i = 0; i < length; i++){
        l2norm += array[i]*array[i];
    }
    l2norm = sqrt(l2norm);
    return l2norm;
}

// Linearly rescale input such that its values are in the range [0, 1].
void linear_conversion(const _myfloat *in, _myfloat *out, int l)
{
    _myfloat a, b;
    _myfloat min = array_min(in, l);
    _myfloat max = array_max(in, l);

    // skip the normalization if max = min
    if( max > min + 0.00001){
        a = (1.0 - 0.0) / (max - min);
        b = -a * min;
        for (int i = 0; i < l; i++)
            out[i] = a * in[i] + b;
    }
    else{
        for (int i = 0; i < l; i++)
            out[i] = in[i];
    }
}

// Compute the x modulus y
_myfloat modulus(_myfloat x, _myfloat y)
{
    _myfloat z = x;
    int n;
    if(z < 0){
        n = (int)((-z)/y)+1;
        z += n*y;
    }
    n = (int)(z/y);
    z -= n*y;
    return z;
}

// Multiply the rotation matric R_alpha with [x,y]^T
void apply_rotation(_myfloat x, _myfloat y, _myfloat *rx, _myfloat *ry, _myfloat alpha)
{
    _myfloat c = cos(alpha);
    _myfloat s = sin(alpha);
    _myfloat tx = c * x - s * y;
    _myfloat ty = s * x + c * y;
    *rx = tx;
    *ry = ty;
}





// EXTRA
void array_mean_and_std(const _myfloat* array, int length, _myfloat* mean, _myfloat* std)
{
    // mean value
    _myfloat sum = 0;
    for(int i = 0; i < length; i++){
        sum += array[i];
    }
    *mean = sum/(_myfloat)length;

    // standard deviation
    sum = 0;
    _myfloat tmp;
    for(int i = 0; i < length; i++){
        tmp = array[i]-(*mean);
        sum += tmp*tmp; 
    }
    *std = sqrt(sum/(_myfloat)length);
}
