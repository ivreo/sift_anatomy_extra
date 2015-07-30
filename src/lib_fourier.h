

#ifndef _LIB_FOURIER_H_
#define _LIB_FOURIER_H_

#include <fftw3.h>
#include "lib_util.h" // for the definition of _myfloat


void laplacian_with_dct(_myfloat* x, _myfloat* y, int w, int h);

void add_gaussian_blur_dct(_myfloat* in, _myfloat* out, int w, int h, _myfloat sigma);


void direct_dct(_myfloat* im, _myfloat* dct, int w, int h);
void inverse_dct(_myfloat* dct, _myfloat* im, int w, int h);


void crop_dct(_myfloat* dctin, _myfloat* dctout, int wi, int hi, int wo, int ho);
void zeropad_dct(_myfloat* dctin, _myfloat* dctout, int wi, int hi, int wo, int ho);
void resample_dct(_myfloat* dctin, _myfloat* dctout, int wi, int hi, int wo, int ho);

void zoom_by_zeropadding_dct(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho);
void subsampling_crop_dct(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho);

void gaussian_on_dct(const _myfloat* dctin, _myfloat* dctout, int w, int h, _myfloat sigma);




/*  DFT  */

void laplacian_with_dft(_myfloat* in, _myfloat* out, int w, int h);

void zoom_by_zeropadding_dft(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho);
void subsampling_crop_dft(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho);

void oversample_DCT(const _myfloat* in, int wi, int hi, _myfloat* out, int wo, int ho);


#endif /* _LIB_FOURIER_H_ */
