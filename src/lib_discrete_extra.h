#ifndef _LIB_DISCRETE_EXTRA_H_
#define _LIB_DISCRETE_EXTRA_H_

#include "lib_util.h"

void compute_laplacian(_myfloat* in, _myfloat* out, int w, int h);
void compute_laplacian_scheme(_myfloat* in, _myfloat* out, int w, int h, _myfloat l);
void subsample_by_intfactor(_myfloat* in, _myfloat* out, int wi, int hi, int factor);


#endif // _LIB_DISCRETE_EXTRA_H_
