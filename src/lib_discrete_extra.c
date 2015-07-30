/*
IPOL SIFT
Copyright (C) 2013, Ives Rey Otero, CMLA ENS Cachan 
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )



#include "lib_discrete_extra.h"

void compute_laplacian(_myfloat* in, _myfloat* out, int w, int h)
{
    for(int i=1; i<h-1; i++){
        for(int j=1; j<w-1; j++){
            out[i*w+j] = in[(i+1)*w+j] + in[(i-1)*w+j] \
                       + in[i*w+j+1] + in[i*w+j-1] \
                       - 4*in[i*w+j];
        }
    }
}


void compute_laplacian_scheme(_myfloat* in, _myfloat* out, int w, int h, _myfloat l)
{
    assert((0<=l)&&(l <=1));
    for(int i=1; i<h-1; i++){
        for(int j=1; j<w-1; j++){
            // \Delta^{+}
            _myfloat lapl_p = in[(i+1)*w+j] + in[(i-1)*w+j] \
                          + in[i*w+j+1] + in[i*w+j-1] \
                          - 4*in[i*w+j];
            // \Delta^{\times}
            _myfloat lapl_t = 0.5*( in[(i+1)*w+(j+1)] + in[(i-1)*w+(j+1)] \
                                + in[(i+1)*w+(j-1)] + in[(i-1)*w+(j-1)] \
                                - 4* in[i*w+j]);
            out[i*w+j] = l*lapl_p + (1-l)*lapl_t;
        }
    }
}



/** @brief sub sampling by an integer factor, keeping sample (0,0) */
void subsample_by_intfactor(_myfloat* in, _myfloat* out, int wi, int hi, int factor)
{
    int wo = wi/factor;
    int ho= hi/factor;
    for(int i=0;i<ho;i++){
        int i_p=factor*i;
        for(int j=0;j<wo;j++){
            int j_p=factor*j;
            out[i*wo+j] = in[i_p*wi+j_p];
        }
    }
}




/** @brief Compute the determinant of the 2D Hessian
 *
 * Hessian schemes : 
 *
 * d^2 u / dx^2 = u(x+1) + u(x-1) -2*u(x)
 *
 * d^2 u / dx dy = ( u(x-1,y-1) + u(x+1,y+1) - u(x-1,y+1) - u(x+1,y-1) )/4.0
 *
 * Det : Hxx*Hyy - Hxy*Hxy
 *
 */
void compute_hessian_determinant(_myfloat* in, _myfloat* out, int w, int h)
{
    for(int i=1; i<h-1; i++){
        for(int j=1; j<w-1; j++){


            _myfloat dxx = in[(i-1)*w+j]
                      -2*in[i*w+j]
                        +in[(i+1)*w+j];
  
            _myfloat dyy = in[i*w+j-1]
                      -2*in[i*w+j]
                        +in[i*w+j+1];

            _myfloat dxy = ( in[(i-1)*w + (j-1)]
                          -in[(i-1)*w + (j+1)]
                          -in[(i+1)*w + (j-1)]
                          +in[(i+1)*w + (j+1)])/4.0;

            out[i*w+j] = dxx*dyy - dxy*dxy;
        }
    }
}




