#ifndef _BSPLINES_H_
#define _BSPLINES_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>




double* image_bspline_decomposition(double *in, int order, int width, int height);
double* precompute_bspline_binomials(int order);
double  image_bspline_reconstruction(double x, double y, int order, double* binomials, double* bspline_decomp_image, int width, int height);





#endif