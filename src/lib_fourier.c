
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <fftw3.h>




#include "lib_fourier.h"
#include "lib_util.h"




static void symmetrize_image(const _myfloat* in, _myfloat* out, int w, int h)
{
    for(int i = 0; i < h; i++){
        for(int j = 0; j < w; j++){
            _myfloat v = in[i*w+j];
            out[i*(2*w)+j] = v;
            out[(2*h-1-i)*(2*w)+j] = v;
            out[i*(2*w)+(2*w-1-j)] = v;
            out[(2*h-1-i)*(2*w)+(2*w-1-j)] = v;
        }
    }
}

static void take_topleftquarter(_myfloat* in, _myfloat* out, int wi, int hi)
{
    for(int i = 0; i < hi/2; i++){
        for(int j = 0; j < wi/2; j++){
            out[i*(wi/2)+j] = in[i*wi+j];
        }
    }
}





//////////////////
void direct_dct(_myfloat* im, _myfloat* dct, int w, int h)
{
#ifdef QUAD
    fftwq_plan r2r;
    r2r = fftwq_plan_r2r_2d(h, w, im,  dct, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    fftwq_execute(r2r);
    _myfloat norm = (_myfloat)(2*w*2*h);
    for(int i=0; i<w*h; i++)
        dct[i] /=norm;
    fftwq_destroy_plan(r2r);
#else
    #ifdef DOUBLE
    fftw_plan r2r;
    //r2r = fftw_plan_r2r_2d(h, w, (double*)im, (double*)dct, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    r2r = fftw_plan_r2r_2d(h, w, im, dct, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(r2r);
    _myfloat norm = (_myfloat)(2*w*2*h);
    for(int i=0; i<w*h; i++)
        dct[i] /=norm;
    fftw_destroy_plan(r2r);
    #else
    fftwf_plan r2r;
    //r2r = fftwf_plan_r2r_2d(h, w, (float*)im, (float*)dct, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    r2r = fftwf_plan_r2r_2d(h, w, im, dct, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
    fftwf_execute(r2r);
    _myfloat norm = (_myfloat)(2*w*2*h);
    for(int i=0; i<w*h; i++)
        dct[i] /=norm;
    fftwf_destroy_plan(r2r);
    #endif
#endif






}


/////////////////////
void inverse_dct(_myfloat* dct, _myfloat* im, int w, int h)
{
#ifdef QUAD
    fftwq_plan r2r;
	r2r = fftwq_plan_r2r_2d(h, w, dct, im, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    fftwq_execute(r2r);
    fftwq_destroy_plan(r2r);
#else
    #ifdef DOUBLE
    fftw_plan r2r;
	r2r = fftw_plan_r2r_2d(h, w, dct, im, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(r2r);
    fftw_destroy_plan(r2r);
    #else
    fftwf_plan r2r;
	r2r = fftwf_plan_r2r_2d(h, w, dct, im, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    fftwf_execute(r2r);
    fftwf_destroy_plan(r2r);
    #endif
#endif

}






/* WARNING :
 *  
 *  DCT and DFT of symmetrized signal are related, however
 *
 *  the complex factor alpha_n^N =  exp(+i\pi n / 2N) such that
 *
 *  DFT(sym(u))_n = alpha_n^N DCT(u)_n
 *
 *  make it impossible to make a zero padding without some sort of translation  
 *
 * The DCT interpolation for the original and zero padded signal are
 *
 *  u(x) = dct_0 + \sum_{n = 1\ldots N} dct_n cos(\pi n /N ( x+1/2))
 *
 *  v(x) = dct_0 + \sum_{n = 1\ldots M} dct_n cos(\pi n /M ( x+1/2))
 *
 *  so we have u(x) = v(y) for
 *
 *     x = N/M y + (N/(2M)) -1/2
 *
 *
 * ****************
 *  For Fourier zoom outs I recommend to zero pad the DFT of summetrized signal instead.
 * ****************
 *
 */
void gaussian_on_dct(const _myfloat* dctin, _myfloat* dctout, int w, int h, _myfloat sigma)
{
    _myfloat s = sigma*sigma*M_PI*M_PI/2.0;
    for(int m=0; m<h; m++){
        _myfloat dm = (_myfloat)m/(_myfloat)h;
        for(int n = 0; n < w; n++){
            _myfloat dn = (_myfloat)n/(_myfloat)w;
            dctout[m*w+n] = dctin[m*w+n]*exp(-s*(dm*dm+dn*dn));
        }
    }
}




void add_gaussian_blur_dct(_myfloat* in, _myfloat* out, int w, int h, _myfloat sigma)
{
    _myfloat* dctin = malloc(w*h*sizeof(_myfloat));
    _myfloat* dctout = malloc(w*h*sizeof(_myfloat));
    direct_dct(in, dctin, w, h);
    gaussian_on_dct(dctin, dctout, w, h, sigma);
    inverse_dct(dctout, out, w, h);
    free(dctin);
    free(dctout);
}











// DECALAGE 
void zeropad_dct(_myfloat* dctin, _myfloat* dctout, int wi, int hi, int wo, int ho)
{
    for(int m = 0; m < ho*wo; m++){
        dctout[m]=0.;
    }
    for(int m = 0; m < hi; m++){
        for(int n = 0; n < wi; n++){
            dctout[m*wo+n] = dctin[m*wi+n];
        }
    }
}




void crop_dct(_myfloat* dctin, _myfloat* dctout, int wi, int hi, int wo, int ho)
{
    assert(ho<=hi);
    assert(wo<=wi);

    for(int m=0; m<ho; m++){
        for(int n=0; n<wo; n++){
            dctout[m*wo+n] = dctin[m*wi+n];
        }
    }
}


void resample_dct(_myfloat* dctin, _myfloat* dctout, int wi, int hi, int wo, int ho)
{
    if((wo>wi)&&(ho>hi)){
        zeropad_dct(dctin, dctout, wi, hi, wo, ho);
    }
    
    if((wo<wi)&&(ho<hi)){
        crop_dct(dctin, dctout, wi, hi, wo, ho);
    }
    
    if((wo==wi)&&(ho==hi)){
        for(int i=0; i<wo*ho; i++){
            dctout[i] = dctin[i];
        }
    }
}





void zoom_by_zeropadding_dct(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho)
{
    _myfloat* dctin  = malloc(wi*hi*sizeof(_myfloat));
    _myfloat* dctout = malloc(wo*ho*sizeof(_myfloat));

    direct_dct(in, dctin, wi, hi);
    zeropad_dct(dctin, dctout, wi, hi, wo, ho);
    inverse_dct(dctout, out, wo, ho);

    free(dctin);
    free(dctout);
}


void subsampling_crop_dct(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho)
{
    _myfloat* dctin  = malloc(wi*hi*sizeof(_myfloat));
    _myfloat* dctout = malloc(wo*ho*sizeof(_myfloat));

    direct_dct(in, dctin, wi, hi);
    crop_dct(dctin, dctout, wi, hi, wo, ho);
    inverse_dct(dctout, out, wo, ho);

    free(dctin);
    free(dctout);
}




#ifdef QUAD
void direct_dft(_myfloat* im, fftwq_complex* dft, int w, int h)
{
    fftwq_plan r2c;
	r2c = fftwq_plan_dft_r2c_2d(h, w, im, dft, FFTW_ESTIMATE);
	fftwq_execute(r2c);
    fftwq_destroy_plan(r2c);
    _myfloat norm = (_myfloat)(w*h);
    for(int i=0; i<(w/2+1)*h; i++){
        dft[i][0] /=norm;
        dft[i][1] /=norm;
    }
}
#else
    #ifdef DOUBLE
void direct_dft(_myfloat* im, fftw_complex* dft, int w, int h)
{
    fftw_plan r2c;
	r2c = fftw_plan_dft_r2c_2d(h, w, im, dft, FFTW_ESTIMATE);
	fftw_execute(r2c);
    fftw_destroy_plan(r2c);
    _myfloat norm = (_myfloat)(w*h);
    for(int i=0; i<(w/2+1)*h; i++){
        dft[i][0] /=norm;
        dft[i][1] /=norm;
    }
}
    #else
void direct_dft(_myfloat* im, fftwf_complex* dft, int w, int h)
{
    fftwf_plan r2c;
	r2c = fftwf_plan_dft_r2c_2d(h, w, im, dft, FFTW_ESTIMATE);
	fftwf_execute(r2c);
    fftwf_destroy_plan(r2c);
    _myfloat norm = (_myfloat)(w*h);
    for(int i=0; i<(w/2+1)*h; i++){
        dft[i][0] /=norm;
        dft[i][1] /=norm;
    }
}
    #endif
#endif



/** 
 *  dft : size (w/2+1)*h 
 */
#ifdef QUAD
void inverse_dft(fftwq_complex* dft, _myfloat* im, int w, int h)
{
    fftwq_plan c2r;
	c2r = fftwq_plan_dft_c2r_2d(h, w, dft, im, FFTW_ESTIMATE);
	fftwq_execute(c2r);
    fftwq_destroy_plan(c2r);
}
#else
    #ifdef DOUBLE
void inverse_dft(fftw_complex* dft, _myfloat* im, int w, int h)
{
    fftw_plan c2r;
	c2r = fftw_plan_dft_c2r_2d(h, w, dft, im, FFTW_ESTIMATE);
	fftw_execute(c2r);
    fftw_destroy_plan(c2r);
}
    #else
void inverse_dft(fftwf_complex* dft, _myfloat* im, int w, int h)
{
    fftwf_plan c2r;
	c2r = fftwf_plan_dft_c2r_2d(h, w, dft, im, FFTW_ESTIMATE);
	fftwf_execute(c2r);
    fftwf_destroy_plan(c2r);
}
    #endif
#endif


// DEGOOOOOLASSE 
#ifdef QUAD
void gaussian_on_dft(fftwq_complex* dftin, fftwq_complex* dftout, int w, int h, _myfloat sigma)
#else
    #ifdef DOUBLE
void gaussian_on_dft(fftw_complex* dftin, fftw_complex* dftout, int w, int h, _myfloat sigma)
    #else
void gaussian_on_dft(fftwf_complex* dftin, fftwf_complex* dftout, int w, int h, _myfloat sigma)
    #endif
#endif
{
    _myfloat s = 2*sigma*sigma*M_PI*M_PI;
	for(int m = 0; m < h; m++){
        int mp = m;
		if (m > -h/2+h-1){mp = m-h;}
        _myfloat dmp = (_myfloat)mp/(_myfloat)h;
		for(int n = 0; n < w/2+1; n++){
			int np = n;
			if (n > -w/2+w-1){ np = n-w;}
            _myfloat dnp = (_myfloat)np/(_myfloat)w;
			_myfloat weight = exp(-s*(dmp*dmp + dnp*dnp));
			dftout[m*(w/2+1)+n][0] = dftin[m*(w/2+1)+n][0] * weight;
			dftout[m*(w/2+1)+n][1] = dftin[m*(w/2+1)+n][1] * weight;
		}
	}
}



    
#ifdef QUAD
void laplacian_on_dft(fftwq_complex* dftin, fftwq_complex* dftout, int w, int h)
#else
    #ifdef DOUBLE
void laplacian_on_dft(fftw_complex* dftin, fftw_complex* dftout, int w, int h)
    #else
void laplacian_on_dft(fftwf_complex* dftin, fftwf_complex* dftout, int w, int h)
    #endif
#endif
{
	for(int m = 0; m < h; m++){
        int mp = m;
		if (m > -h/2+h-1){mp = m-h;}
        _myfloat dmp = (_myfloat)mp/(_myfloat)h;
		for(int n = 0; n < w/2+1; n++){
			int np = n;
			if (n > -w/2+w-1){ np = n-w;}
            _myfloat dnp = (_myfloat)np/(_myfloat)w;
			_myfloat weight = 4*M_PI*M_PI*(dmp*dmp + dnp*dnp);
			dftout[m*(w/2+1)+n][0] = dftin[m*(w/2+1)+n][0] * weight;
			dftout[m*(w/2+1)+n][1] = dftin[m*(w/2+1)+n][1] * weight;
		}
	}
}

#ifdef QUAD
void laplacian_with_dft(_myfloat* in, _myfloat* out, int w, int h)
{

	fftwq_complex* dftin;
	fftwq_complex* dftout;
    dftin  = (fftwq_complex*)fftwq_malloc((w/2+1)*h*sizeof(fftwq_complex));
	dftout = (fftwq_complex*)fftwq_malloc((w/2+1)*h*sizeof(fftwq_complex));
    
    direct_dft(in, dftin, w, h);
    laplacian_on_dft(dftin, dftout, w, h);
    inverse_dft(dftout, out, w, h);
	
    fftwq_free(dftin);
	fftwq_free(dftout);
}
#else
    #ifdef DOUBLE
void laplacian_with_dft(_myfloat* in, _myfloat* out, int w, int h)
{

	fftw_complex* dftin;
	fftw_complex* dftout;
    dftin  = (fftw_complex*)fftw_malloc((w/2+1)*h*sizeof(fftw_complex));
	dftout = (fftw_complex*)fftw_malloc((w/2+1)*h*sizeof(fftw_complex));
    
    direct_dft(in, dftin, w, h);
    laplacian_on_dft(dftin, dftout, w, h);
    inverse_dft(dftout, out, w, h);
	
    fftw_free(dftin);
	fftw_free(dftout);
}
    #else 
void laplacian_with_dft(_myfloat* in, _myfloat* out, int w, int h)
{

	fftwf_complex* dftin;
	fftwf_complex* dftout;
    dftin  = (fftwf_complex*)fftwf_malloc((w/2+1)*h*sizeof(fftwf_complex));
	dftout = (fftwf_complex*)fftwf_malloc((w/2+1)*h*sizeof(fftwf_complex));
    
    direct_dft(in, dftin, w, h);
    laplacian_on_dft(dftin, dftout, w, h);
    inverse_dft(dftout, out, w, h);
	
    fftwf_free(dftin);
	fftwf_free(dftout);
}
    #endif
#endif






/* Using DFT and symmetrization */
void laplacian_with_dct(_myfloat* x, _myfloat* y, int w, int h)
{
    _myfloat* xx = malloc(2*w*2*h*sizeof(*xx));
    _myfloat* yy = malloc(2*w*2*h*sizeof(*yy));
    symmetrize_image(x, xx, w, h);
    laplacian_with_dft(xx, yy, 2*w, 2*h);
    take_topleftquarter(yy, y, 2*w, 2*h);
    free(xx);
    free(yy);
}




#ifdef QUAD
void zeropad_dft_bis(fftwq_complex* dftin, fftwq_complex* dftout, int wi, int hi, int wo, int ho)
#else
#ifdef DOUBLE
void zeropad_dft_bis(fftw_complex* dftin, fftw_complex* dftout, int wi, int hi, int wo, int ho)
#else
void zeropad_dft_bis(fftwf_complex* dftin, fftwf_complex* dftout, int wi, int hi, int wo, int ho)
#endif
#endif
{
    assert(ho>=hi);
    assert(wo>=wi);

    for(int i = 0; i < (wo/2+1)*ho; i++){
        dftout[i][0] = 0.;
        dftout[i][1] = 0.;
    }


    for(int i = 0; i < hi; i++){
        int a = i;
        if (a > -hi/2+(hi-1)){
            a = a - hi;
        }
        for(int j = 0; j < wi/2; j++){
            int b = j;
            if (b > -wi/2+(wi-1)){b = b - wi;}

            int c = a;
            int d = b;
            if(c<0){c = c+ho;}
            if(d<0){d = c+wo;}

            dftout[c*(wo/2+1)+d][0] = dftin[i*(wi/2+1)+j][0];
            dftout[c*(wo/2+1)+d][1] = dftin[i*(wi/2+1)+j][1];

        }
    }
}


///////////////////////////////
#ifdef QUAD
void zoom_by_zeropadding_dft(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho)
{
	fftwq_complex* dftin;
	fftwq_complex* dftout;
    dftin  = (fftwq_complex*)fftwf_malloc((wi/2+1)*hi*sizeof(fftwq_complex));
	dftout = (fftwq_complex*)fftwf_malloc((wo/2+1)*ho*sizeof(fftwq_complex));
    
    direct_dft(in, dftin, wi, hi);
    zeropad_dft_bis(dftin, dftout, wi, hi, wo, ho);
    inverse_dft(dftout, out, wo, ho);

    fftwq_free(dftin);
    fftwq_free(dftout);
}
#else
    #ifdef DOUBLE
void zoom_by_zeropadding_dft(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho)
{
	fftw_complex* dftin;
	fftw_complex* dftout;
    dftin  = (fftw_complex*)fftw_malloc((wi/2+1)*hi*sizeof(fftw_complex));
	dftout = (fftw_complex*)fftw_malloc((wo/2+1)*ho*sizeof(fftw_complex));
    
    direct_dft(in, dftin, wi, hi);
    zeropad_dft_bis(dftin, dftout, wi, hi, wo, ho);
    inverse_dft(dftout, out, wo, ho);

    fftw_free(dftin);
    fftw_free(dftout);
}
    #else
void zoom_by_zeropadding_dft(_myfloat* in, _myfloat* out, int wi, int hi, int wo, int ho)
{
	fftwf_complex* dftin;
	fftwf_complex* dftout;
    dftin  = (fftwf_complex*)fftwf_malloc((wi/2+1)*hi*sizeof(fftwf_complex));
	dftout = (fftwf_complex*)fftwf_malloc((wo/2+1)*ho*sizeof(fftwf_complex));
    
    direct_dft(in, dftin, wi, hi);
    zeropad_dft_bis(dftin, dftout, wi, hi, wo, ho);
    inverse_dft(dftout, out, wo, ho);

    fftwf_free(dftin);
    fftwf_free(dftout);
}
    #endif
#endif





// zeropad_DCT mais avec un zeropadding
void oversample_DCT(const _myfloat* in, int wi, int hi, _myfloat* out, int wo, int ho)
{
    _myfloat* inin   = malloc(2*wi*2*hi*sizeof(*inin));
    _myfloat* outout = malloc(2*wo*2*ho*sizeof(*outout));

    symmetrize_image(in, inin, wi, hi);
    zoom_by_zeropadding_dft(inin, outout, 2*wi, 2*hi, 2*wo, 2*ho);
    take_topleftquarter(outout, out, 2*wo,2*ho);
    
    free(inin);
    free(outout);
}

