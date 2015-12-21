#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "lib_description.h"
#include "lib_discrete.h"
#include "lib_sift_anatomy.h"
#include "lib_util.h"

#include "lib_fourier.h"
#include "lib_bsplines.h"
#include "lib_discrete_extra.h"
#include "lib_dense_anatomy.h"

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )
#define ABS(x) ((x)<0?-(x):(x))
#define EPSILON 0

#include "lib_io_scalespace.h"

void oversample_with_flag(const _myfloat* in, int win, int hin, _myfloat* out, int wout, int hout, _myfloat delta_out, int flag_interp);





/** @brief Compute the Gaussian scale-space and the DoG scale-space
 *
 *  - all images are computed from the input image via successive interpolation
 *  (optionally bspline), DCT blur, and subsampling
 *
 *  - for the DoG, each image of the difference is also computed from the
 *  original image (the semi-group is never used).
 *
 *  - note that for the Gaussian scale-space, the supplementary scale won"t be
 *  computed, it is never used anyway since the DoG computed is already
 *  computed, the Gaussian scale-space will only be used
 *
 */
static void scalespace_compute_dense_ss_and_dog(struct sift_scalespace* ss, // Gaussian scale-space
                                     struct sift_scalespace* dd, // DoG
                                     _myfloat* in,
                                     int w_in,
                                     int h_in,
                                     _myfloat sigma_in,
                                     _myfloat k,
                                     int flag_interp)
{

    /* size of the interpolated image */
    _myfloat delta = ss->octaves[0]->delta;
    int w_min = ss->octaves[0]->w;
    int h_min = ss->octaves[0]->h;

    /* checking scale-space definition consistance*/
    assert(w_min == (int)(w_in/delta));
    assert(h_min == (int)(h_in/delta));
    assert(w_min >= w_in);
    assert(h_min >= h_in); // always interpolate

    // Compute the interpolated image (bsplines)
    _myfloat* tmp = xmalloc(w_min*h_min * sizeof(*tmp));
    oversample_with_flag(in, w_in, h_in, tmp, w_min, h_min, delta, flag_interp);
    _myfloat* tmpM = xmalloc(w_min*h_min*sizeof(*tmpM));
    _myfloat* tmpP = xmalloc(w_min*h_min*sizeof(*tmpP));

    // ... and from it, compute the scale-spaces
    for(int o = 0; o < ss->nOct; o++){
        struct octa* d_oct = dd->octaves[o];
        struct octa* oct = ss->octaves[o];
        int w = oct->w;
        int h = oct->h;
        for(int s = 0; s < d_oct->nSca; s++){ /* add blur to previous image in the stack */

            _myfloat sigma = oct->sigmas[s];

            _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta;
            _myfloat sigP = sqrt( k * k * sigma*sigma - sigma_in*sigma_in)/delta;

            add_gaussian_blur_dct(tmp, tmpM, w_min, h_min, sigM);
            add_gaussian_blur_dct(tmp, tmpP, w_min, h_min, sigP);

            _myfloat* imM = xmalloc(w*h*sizeof(*imM));
            _myfloat* imP = xmalloc(w*h*sizeof(*imP));

            int subfactor = pow(2, o);
            subsample_by_intfactor(tmpM, imM, w_min, h_min, subfactor);
            subsample_by_intfactor(tmpP, imP, w_min, h_min, subfactor);

            for(int i = 0; i < w*h; i++){
                // for DoG scale-space
                d_oct->imStack[s*w*h+i] = imP[i] - imM[i];
                // for Gaussian scale-space
                oct->imStack[s*w*h+i] = imM[i];
            }
            xfree(imM);
            xfree(imP);
        }
    }
    xfree(tmp);
    xfree(tmpM);
    xfree(tmpP);
}









// To avoid descriptor computation
void copy_all_keys(struct sift_keypoints* kA, struct sift_keypoints* kB)
{
    for(int i = 0; i < kA->size; i++){
        struct keypoint* onekey = sift_malloc_keypoint_from_model_and_copy(kA->list[i]);
        sift_add_keypoint_to_list(onekey, kB);
    }
}








static void gaussian_blur_with_flag(_myfloat* in, _myfloat* out, int w, int h, _myfloat sigma, int flag_dct);


static void update_and_copy_keypoints_list(struct sift_keypoints *kin,
                                          int o, int s,
                                          struct sift_keypoints *kout)
{
    for(int ik = 0; ik < kin->size; ik++){
        struct keypoint* key = kin->list[ik];
        key->o = o;
        key->s = s;
        struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
        sift_add_keypoint_to_list(copy, kout);
    }
}






struct sift_scalespace* sift_malloc_scalespace_lowe_floatnspo(int nOct,   /* # of octaves            */
                                              //      int nSca,   /* # of scales of detection (excluding the 3 auxiliary scales */
                                                    _myfloat fnspo,
                                                    int im_w, int im_h,  /* # input image dimension */
                                                    _myfloat delta_min,   /* minimal inter distance sample     */
                                                    _myfloat sigma_min)   /* minimal scale in each octave (relatively to the sampling rate) */
{

    int nSca = ceil(fnspo);

    int* ws  = xmalloc(nOct*sizeof(int));
    int* hs  = xmalloc(nOct*sizeof(int));
    int* nScas      = xmalloc(nOct*sizeof(int));
    _myfloat* deltas  = xmalloc(nOct*sizeof(_myfloat));
    _myfloat** sigmas = xmalloc(nOct*sizeof(_myfloat*));    /* nSca+3 = nSca + 2 extra scales to search 3d extrema + 1 extra scale to compute DoG */
    assert(delta_min <=1);
    deltas[0] = delta_min;
    hs[0] = (int)(im_h/delta_min);
    ws[0] = (int)(im_w/delta_min);
    for(int o=1;o<nOct;o++){
        ws[o] = ws[o-1]/2;     /*integer division*/
        hs[o] = hs[o-1]/2;
        deltas[o] = delta_min * pow(2,o); // deltas[o-1]*2.0;
    }

    
    for(int o=0;o<nOct;o++){
        nScas[o] = nSca+3;  /* 3 extra images in the stack, 1 for dog computation and 2 for 3d discrete extrema definition*/
        sigmas[o] = xmalloc(nScas[o]*sizeof(_myfloat));
        for(int s=0;s<nSca+3;s++){ /* nSca images + 3 auxiliary images*/
            sigmas[o][s] = sigma_min*pow(2.0, o + (_myfloat)s/fnspo);
        }
    }
    struct sift_scalespace* scalespace = sift_malloc_scalespace(nOct, deltas, ws, hs, nScas, sigmas);   
    xfree(deltas);
    xfree(ws);
    xfree(hs);
    xfree(nScas);
    for(int o=0;o<nOct;o++)
        xfree(sigmas[o]);
    xfree(sigmas);    

    return scalespace;
}



// A copy of convert_threshold (but uses dog_nspo instead of n_spo)
_myfloat convert_threshold_DoG(const struct sift_parameters* p)
{
    // converting the threshold to make it consistent
    _myfloat k_nspo = pow(2, (_myfloat)(p->dog_nspo));
    _myfloat k_3 =  pow(2, 3.0);

    _myfloat thresh = (k_nspo - 1) / (k_3 - 1) * p->C_DoG ;
    return thresh;
}



struct sift_keypoints* sift_anatomy_dense(_myfloat* x, int w, int h,
                                          struct sift_parameters* p,
                                          struct sift_scalespace* ss[4],
                                          struct sift_keypoints* kk[6],
                                          int flag_interp)
{
    struct sift_keypoints* k = sift_malloc_keypoints();
    int n_oct = number_of_octaves(w, h, p);

    // normalizing the threshold
    _myfloat thresh = convert_threshold_DoG(p);

    /** MEMORY ALLOCATION **/
    /** scale-space structure */
    struct sift_scalespace* s   = sift_malloc_scalespace_lowe_floatnspo(n_oct, p->fnspo, w, h, p->delta_min, p->sigma_min);
    struct sift_scalespace* d   = sift_malloc_scalespace_dog_from_scalespace(s);
    /** list-of-keypoints (Already allocated) */
    struct sift_keypoints* kA   = kk[0];  /* 3D (discrete) extrema   */
    struct sift_keypoints* kB   = kk[1];  /* passing the threshold on DoG  */
    struct sift_keypoints* kC   = kk[2];  /* interpolated 3D extrema (continuous) */
    struct sift_keypoints* kD   = kk[3];  /* passing the threshold on DoG  */
    struct sift_keypoints* kE   = kk[4];  /* passing OnEdge filter */
    struct sift_keypoints* kF   = kk[5];  /* image border */

    /** SCALESPACES COMPUTATION **********************************************/
    _myfloat kfact = pow(2, 1.0/(_myfloat)(p->dog_nspo)); // equivalent definition of the DoG operator
    scalespace_compute_dense_ss_and_dog(s, d, x, w, h, p->sigma_in, kfact, flag_interp);

    
    /** KEYPOINT DETECTION ***************************************************/
    keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
    keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
    keypoints_interpolate_position_controlled(d, kB, kC, p->itermax, p->ofstMax_X, p->ofstMax_S, p->flag_jumpinscale, p->fnspo, p->sigma_min);
    keypoints_discard_with_low_response(kC, kD, thresh);
    keypoints_compute_edge_response(d,kD);
    keypoints_discard_on_edge( kD, kE, (p->C_edge+1)*(p->C_edge+1)/p->C_edge );
    keypoints_discard_near_the_border(kE, kF, w, h, 1.0);  // should have a parameter


    /** KEYPOINT DESCRIPTION ***************************************************/
    struct sift_scalespace* sx  = sift_malloc_scalespace_from_model(s);
    struct sift_scalespace* sy  = sift_malloc_scalespace_from_model(s);
    scalespace_compute_gradient(s,sx,sy); /* Pre-computes gradient scale-space */
    keypoints_attribute_orientations(sx, sy, kF, k, p->n_bins,p->lambda_ori,p->t);
    keypoints_attribute_descriptors(sx, sy, k, p->n_hist,p->n_ori,p->lambda_descr);

    /** scalespace structures*/
    ss[0] = s;
    ss[1] = d;
    ss[2] = sx;
    ss[3] = sy;

    return k;
}




// Minimal memory allocation for the gradual computation of one octave (Computing 3 scales)
struct sift_scalespace* sift_malloc_scalespace_dog_gradual_floatnspo(int ObjnOct,
                                                       _myfloat ObjnSca,
                                                       int im_w, int im_h,
                                                       _myfloat delta_min, _myfloat sigma_min,
                                                       int curr_o, int curr_s) // current coordinates
{
    // only compute a octave of three consecutive scales.
    int nOct = 1;
    int nSca = 3; //pour le dog
    // normal
    int* ws = xmalloc(nOct*sizeof(int));
    int* hs = xmalloc(nOct*sizeof(int));
    int* nScas = xmalloc(nOct*sizeof(int));
    _myfloat* deltas = xmalloc(nOct*sizeof(_myfloat));
    _myfloat** sigmas = xmalloc(nOct*sizeof(_myfloat*));

    deltas[0] = delta_min * pow(2, curr_o);
    ws[0] = (int)(im_w / (delta_min * pow(2, curr_o)) );
    hs[0] = (int)(im_h / (delta_min * pow(2, curr_o)) );
    _myfloat fsigma_min = sigma_min * pow(2, curr_o);
    nScas[0] = nSca;
    sigmas[0] = xmalloc(nScas[0]*sizeof(_myfloat));
    for(int s = -1; s <= 1; s++){
        sigmas[0][s+1] = fsigma_min*pow(2.0,(_myfloat)(s+curr_s)/ObjnSca);
    }
    struct sift_scalespace* scalespace = sift_malloc_scalespace(nOct, deltas, ws, hs, nScas, sigmas);
    xfree(deltas);
    xfree(ws);
    xfree(hs);
    xfree(nScas);
    xfree(sigmas[0]);
    xfree(sigmas);
    return scalespace;
}




/** @brief Compute Gaussian convolution 
 *    flag_dct values:
 *      0 : discrete convolution
 *      1 : continuous convolution of the DCT interpolate
 *
 */
static void gaussian_blur_with_flag(_myfloat* in, _myfloat* out, int w, int h, _myfloat sigma, int flag_dct)
{
    if (flag_dct == 1){
        add_gaussian_blur_dct(in, out, w, h, sigma);
    }
    else{
        sift_add_gaussian_blur(in, out, w, h, sigma);
    }
}




/** @brief Oversample image with various interpolation methods
 *    flag_interp values:
 *      0 : bilinear
 *      1 : DCT
 *      3,5,7,9,11 : Bsplines
 *
 */
void oversample_with_flag(const _myfloat* in, int win, int hin,
                          _myfloat* out, int wout, int hout,
                          _myfloat delta_out, int flag_interp)
{
    if (flag_interp == 1){
        oversample_DCT(in, win, hin, out, wout, hout);
    }
    else if ((flag_interp >=2) && (flag_interp <= 11)){
        int order = flag_interp;
        oversample_bsplines(in, win, hin,
                           out, wout, hout,
                           delta_out, order);
    }
    else{
        sift_oversample_bilin(in, win, hin, out, wout, hout, delta_out);
    }
}



/** @brief compute the interpolated image used to compute the scale-space
 *
 *
 *
 */
_myfloat* compute_interpolated(const _myfloat* in, int w, int h,
                               struct sift_parameters* p,
                               int* w_seed, int* h_seed,
                               int flag_interp)
{
    _myfloat delta = p->delta_min;
    int wout = (int)(w/delta);
    int hout = (int)(h/delta);
    _myfloat* out = xmalloc(wout*hout*sizeof(*out));
    // oversampling
    oversample_with_flag(in, w, h, out, wout, hout, delta, flag_interp);
    *w_seed = wout;
    *h_seed = hout;
    return(out);
}





// Computes a one octave scale-space
//   - images are computed from the interpolated image
void compute_controlled_nosemi_dog_from_interpolated(_myfloat* in, int w_in, int h_in,
                                                     _myfloat delta_in, _myfloat sigma_in,
                                                     struct sift_scalespace *d,
                                                     _myfloat k,
                                                     int flag_dct,
                                                     int subfactor)
{
    _myfloat* tmp = xmalloc(w_in*h_in*sizeof(*tmp));
    struct octa* d_oct = d->octaves[0];
    int ns = d_oct->nSca;
    int w = d_oct->w;
    int h = d_oct->h;
    for(int s = 0; s < ns; s++){
        // image lower scale (in DoG)
        _myfloat sigma = d_oct->sigmas[s] ;
        _myfloat* imM = xmalloc(w*h*sizeof(*imM));
        _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta_in;
        gaussian_blur_with_flag(in, tmp, w_in, h_in, sigM, flag_dct);
        subsample_by_intfactor(tmp, imM, w_in, h_in, subfactor);
        // image larger scale (in DoG)
        _myfloat* imP = xmalloc(w*h*sizeof(*imP));
        _myfloat sigP = sqrt( k * k * sigma*sigma - sigma_in*sigma_in)/delta_in;
        gaussian_blur_with_flag(in, tmp, w_in, h_in, sigP, flag_dct);
        subsample_by_intfactor(tmp, imP, w_in, h_in, subfactor);
        // difference
        for(int i = 0; i < w*h; i++){
            d_oct->imStack[s*w*h+i] = imP[i] - imM[i];
        }
        xfree(imP);
        xfree(imM);
    }
    xfree(tmp);
}





void update_controlled_nosemi_dog_from_interpolated(_myfloat* in, int w_in, int h_in,
                                                    _myfloat delta_in, _myfloat sigma_in,
                                                    struct sift_scalespace *d,
                                                    _myfloat k,
                                                    int flag_dct,
                                                    int subfactor)
{
    _myfloat* tmp = xmalloc(w_in*h_in*sizeof(*tmp));
    struct octa* d_oct = d->octaves[0];
    int w = d_oct->w;
    int h = d_oct->h;

    // COPY previously computed images
    for(int i = 0; i < 2*w*h; i++){
        d_oct->imStack[i] = d_oct->imStack[i+w*h];
    }

    // ADD the new image
    // image lower scale (in DoG)
    _myfloat new_sigma = d_oct->sigmas[2];  // supposing it has been updated before
    _myfloat* imM = xmalloc(w*h*sizeof(*imM));
    _myfloat sigM = sqrt(new_sigma*new_sigma - sigma_in*sigma_in)/delta_in;
    gaussian_blur_with_flag(in, tmp, w_in, h_in, sigM, flag_dct);
    subsample_by_intfactor(tmp, imM, w_in, h_in, subfactor);

    // image higher scale (in DoG)
    _myfloat* imP = xmalloc(w*h*sizeof(*imP));
    _myfloat sigP = sqrt( k * k * new_sigma*new_sigma - sigma_in*sigma_in)/delta_in;
    gaussian_blur_with_flag(in, tmp, w_in, h_in, sigP, flag_dct);
    subsample_by_intfactor(tmp, imP, w_in, h_in, subfactor);

    // difference of gaussians
    for(int i = 0; i < w*h; i++){
        d_oct->imStack[2*w*h+i] = imP[i] - imM[i];
    }

    xfree(imP);
    xfree(imM);
    xfree(tmp);
}



/** @brief Computes one octave scale-spaces (DoG and Gaussian)
 *   all images are computed from the interpolated image
 *
 */
static void compute_dog_and_scalespace_from_interpolated(_myfloat* in, int w_in, int h_in,
                                                  _myfloat delta_in, _myfloat sigma_in,
                                                  struct sift_scalespace *d,
                                                  struct sift_scalespace *ss,
                                                  _myfloat k,
                                                  int flag_dct,
                                                  int subfactor)
{
    _myfloat* tmp = xmalloc(w_in*h_in*sizeof(*tmp));
    struct octa* d_oct = d->octaves[0];
    struct octa* ss_oct = ss->octaves[0];
    int ns = d_oct->nSca;
    int w = d_oct->w;
    int h = d_oct->h;
    for(int s = 0; s < ns; s++){
        // image lower scale (in DoG)
        _myfloat sigma = d_oct->sigmas[s] ;
        _myfloat* imM = xmalloc(w*h*sizeof(*imM));
        _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta_in;
        gaussian_blur_with_flag(in, tmp, w_in, h_in, sigM, flag_dct);
        subsample_by_intfactor(tmp, imM, w_in, h_in, subfactor);
        // image larger scale (in DoG)
        _myfloat* imP = xmalloc(w*h*sizeof(*imP));
        _myfloat sigP = sqrt( k * k * sigma*sigma - sigma_in*sigma_in)/delta_in;
        gaussian_blur_with_flag(in, tmp, w_in, h_in, sigP, flag_dct);
        subsample_by_intfactor(tmp, imP, w_in, h_in, subfactor);
        // difference
        for(int i = 0; i < w*h; i++){
            d_oct->imStack[s*w*h+i] = imP[i] - imM[i];
            ss_oct->imStack[s*w*h+i] = imM[i];
        }
        xfree(imP);
        xfree(imM);
    }
    xfree(tmp);
}





static void update_dog_and_scalespace_from_interpolated(_myfloat* in, int w_in, int h_in,
                                                    _myfloat delta_in, _myfloat sigma_in,
                                                    struct sift_scalespace *d,
                                                    struct sift_scalespace *ss,
                                                    _myfloat k,
                                                    int flag_dct,
                                                    int subfactor)
{
    _myfloat* tmp = xmalloc(w_in*h_in*sizeof(*tmp));
    struct octa* d_oct = d->octaves[0];
    struct octa* ss_oct = ss->octaves[0];
    int w = d_oct->w;
    int h = d_oct->h;

    // COPY previously computed images
    for(int i = 0; i < 2*w*h; i++){
        d_oct->imStack[i] = d_oct->imStack[i+w*h];
        ss_oct->imStack[i] = ss_oct->imStack[i+w*h];
    }

    // ADD the new image
    // image lower scale (in DoG)
    _myfloat new_sigma = d_oct->sigmas[2];  // supposing it has been updated before
    _myfloat* imM = xmalloc(w*h*sizeof(*imM));
    _myfloat sigM = sqrt(new_sigma*new_sigma - sigma_in*sigma_in)/delta_in;
    gaussian_blur_with_flag(in, tmp, w_in, h_in, sigM, flag_dct);
    subsample_by_intfactor(tmp, imM, w_in, h_in, subfactor);

    // image higher scale (in DoG)
    _myfloat* imP = xmalloc(w*h*sizeof(*imP));
    _myfloat sigP = sqrt( k * k * new_sigma*new_sigma - sigma_in*sigma_in)/delta_in;
    gaussian_blur_with_flag(in, tmp, w_in, h_in, sigP, flag_dct);
    subsample_by_intfactor(tmp, imP, w_in, h_in, subfactor);

    // difference of gaussians
    for(int i = 0; i < w*h; i++){
        d_oct->imStack[2*w*h+i] = imP[i] - imM[i];
        ss_oct->imStack[2*w*h+i] = imM[i];
    }

    xfree(imP);
    xfree(imM);
    xfree(tmp);
}












/*
 *
 *  First compute the seed image and then compute the scalespace 
 *
 *  Just detection - no thresholds
 *
 */
struct sift_keypoints* sift_anatomy_gradual_old(_myfloat* x, int w, int h,
                                            struct sift_parameters* p,
                                            struct sift_keypoints* kk[6],
                                            int flag_dct,
                                            int flag_interp)
{
    struct sift_keypoints* k = sift_malloc_keypoints();

    int n_oct = number_of_octaves(w, h, p);

    //_myfloat thresh = convert_threshold_DoG(p);

    int w_seed, h_seed;
    _myfloat* seed = compute_interpolated(x, w, h, p, &w_seed, &h_seed, flag_interp);

    for(int o = 0; o < n_oct; o++){

        struct sift_scalespace* d = sift_malloc_scalespace_dog_gradual_floatnspo(n_oct, p->fnspo ,w,h, p->delta_min,p->sigma_min, o, 1);
        int subfactor = pow(2, o);
        int nscales = ceil(p->fnspo);
        for(int is = 0; is < nscales; is++){ // current scale where the keypoints are detected

            /** list-of-keypoints (Already allocated) */
            struct sift_keypoints* kA = sift_malloc_keypoints();
            struct sift_keypoints* ktmp = sift_malloc_keypoints(); // to update the octave and scale indices

            /** KEYPOINT DETECTION ***************************************************/
            _myfloat factDoG = pow(2, 1.0/(_myfloat)(p->dog_nspo)); // the equivalent factor in the definition of the DoG operator
            if(is == 0){
                compute_controlled_nosemi_dog_from_interpolated(seed, w_seed, h_seed, p->delta_min, p->sigma_in,
                                                                    d, factDoG, flag_dct, subfactor);
            }
            else{
                // update scales
                for(int s = 0 ; s < 3 ; s++){
                    d->octaves[0]->sigmas[s] = p->sigma_min*pow(2.0, o + (_myfloat)(is+s) / p->fnspo);
                }
                // the missing image
                update_controlled_nosemi_dog_from_interpolated(seed, w_seed, h_seed, p->delta_min, p->sigma_in,
                                                                   d, factDoG, flag_dct, subfactor);
            }

            if (p->itermax == 0){  // All discrete extrema
                keypoints_find_3d_discrete_extrema_epsilon(d, ktmp, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
            }else{
                keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
                keypoints_interpolate_position_controlled(d, kA, ktmp, p->itermax, p->ofstMax_X, p->ofstMax_S, 0, p->fnspo, pow(2,o)*p->sigma_min);
            }

            // 20150612
            // update the index of octave and scale (by default all keypoints are in scale 1 (gradual implementation)).
            for(int ik = 0; ik < ktmp->size; ik++){
                struct keypoint* key = ktmp->list[ik];
                key->o = o;
                key->s = is+1;
                struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
                sift_add_keypoint_to_list(copy, k);
            }

            sift_free_keypoints(kA);
            sift_free_keypoints(ktmp);
        }
        sift_free_scalespace(d);
    }
    free(seed);
    return k;
}



struct sift_keypoints* sift_anatomy_gradual(_myfloat* x, int w, int h,
                                                  struct sift_parameters* p,
                                                  struct sift_keypoints* kk[6],
                                                  int flag_interp)
{
    struct sift_keypoints* k = sift_malloc_keypoints();
    int n_oct = number_of_octaves(w, h, p);

    // normalizing the threshold
    _myfloat thresh = convert_threshold_DoG(p);

    int w_seed, h_seed;
    _myfloat* seed = compute_interpolated(x, w, h, p, &w_seed, &h_seed, flag_interp);

    // at this stage, one could free the image 

    for(int o = 0; o < n_oct; o++){

        struct sift_scalespace* d = sift_malloc_scalespace_dog_gradual_floatnspo(n_oct, p->fnspo ,w,h, p->delta_min,p->sigma_min, o, 1);
        struct sift_scalespace* ss  = sift_malloc_scalespace_from_model(d);
        int subfactor = pow(2, o);
        int nscales = ceil(p->fnspo);
        for(int is = 0; is < nscales; is++){ // current scale where the keypoints are detected


            /** list-of-keypoints (Already allocated) */
            struct sift_keypoints* kA = sift_malloc_keypoints();  /* 3D (discrete) extrema   */
            struct sift_keypoints* kB = sift_malloc_keypoints();  /* passing the threshold on DoG  */
            struct sift_keypoints* kC = sift_malloc_keypoints();  /* interpolated 3D extrema (continuous) */
            struct sift_keypoints* kD = sift_malloc_keypoints();  /* passing the threshold on DoG  */
            struct sift_keypoints* kE = sift_malloc_keypoints();  /* passing OnEdge filter */
            struct sift_keypoints* kF = sift_malloc_keypoints();  /* image border */
            struct sift_keypoints* kG = sift_malloc_keypoints();  /* described keypoints */
    
            /** SCALESPACES COMPUTATION **********************************************/
            _myfloat factDoG = pow(2, 1.0/(_myfloat)(p->dog_nspo)); 
            if(is == 0){
                // all three images are computed
                compute_dog_and_scalespace_from_interpolated(seed, w_seed, h_seed, p->delta_min, p->sigma_in,
                                                                    d, ss, factDoG, 1, subfactor);
            }
            else{
                // update the scales
                for(int s = 0 ; s < 3 ; s++){
                    d->octaves[0]->sigmas[s] = p->sigma_min*pow(2.0, o + (_myfloat)(is+s) / p->fnspo);
                    ss->octaves[0]->sigmas[s] = p->sigma_min*pow(2.0, o + (_myfloat)(is+s) / p->fnspo);
                }
                // copy the already computed images and compute the missing image
                update_dog_and_scalespace_from_interpolated(seed, w_seed, h_seed, p->delta_min, p->sigma_in,
                                                                   d, ss, factDoG, 1, subfactor);
            }

            /** KEYPOINT DETECTION ***************************************************/
            keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
            keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
            keypoints_interpolate_position_controlled(d, kB, kC, p->itermax, p->ofstMax_X, p->ofstMax_S, 0, p->fnspo, pow(2,o)*p->sigma_min);
            keypoints_discard_with_low_response(kC, kD, thresh);
            keypoints_compute_edge_response(d,kD);
            keypoints_discard_on_edge( kD, kE, (p->C_edge+1)*(p->C_edge+1)/p->C_edge );
            keypoints_discard_near_the_border(kE, kF, w, h, 1.0);


            /** KEYPOINT DESCRIPTION ***************************************************/
            if (1){
                struct sift_scalespace* sx  = sift_malloc_scalespace_from_model(ss);
                struct sift_scalespace* sy  = sift_malloc_scalespace_from_model(ss);
                scalespace_compute_gradient(ss,sx,sy); /* Pre-computes gradient scale-space */
                keypoints_attribute_orientations(sx, sy, kF, kG, p->n_bins,p->lambda_ori,p->t);
                keypoints_attribute_descriptors(sx, sy, kG, p->n_hist,p->n_ori,p->lambda_descr);
                update_and_copy_keypoints_list(kG, o, is+1, k);
                sift_free_scalespace(sx);
                sift_free_scalespace(sy);
            }
            else{
                update_and_copy_keypoints_list(kF, o, is+1, k);
            }

            /** Update the index of octave and scale  */
            update_and_copy_keypoints_list(kA, o, is+1, kk[0]);
            update_and_copy_keypoints_list(kB, o, is+1, kk[1]);
            update_and_copy_keypoints_list(kC, o, is+1, kk[2]);
            update_and_copy_keypoints_list(kD, o, is+1, kk[3]);
            update_and_copy_keypoints_list(kE, o, is+1, kk[4]);
            update_and_copy_keypoints_list(kF, o, is+1, kk[5]);
            sift_free_keypoints(kA);
            sift_free_keypoints(kB);
            sift_free_keypoints(kC);
            sift_free_keypoints(kD);
            sift_free_keypoints(kE);
            sift_free_keypoints(kF);
            sift_free_keypoints(kG);
        }
        sift_free_scalespace(d);
        sift_free_scalespace(ss);
    }
    free(seed);
    return k;
}


//    l'introduire dans le corps de la fonction

// void sift_write_scalespace_binary_file(FILE* fp, struct sift_scalespace* ss)
// {
//     fwrite(&ss->nOct, sizeof(int), 1, fp);
//     for(int i = 0; i < ss->nOct; i++){
//         struct octa* o = ss->octaves[i];
//         fwrite(&o->delta, sizeof(_myfloat), 1, fp);
//         fwrite(&o->w, sizeof(int), 1, fp);
//         fwrite(&o->h, sizeof(int), 1, fp);
//         fwrite(&o->nSca, sizeof(int), 1, fp);
//         fwrite(o->sigmas, sizeof(_myfloat), o->nSca, fp);
//     }
//     for(int i = 0; i < ss->nOct; i++){
//         struct octa* o = ss->octaves[i];
//         fwrite(o->imStack, sizeof(_myfloat), o->nSca*o->h*o->w, fp);
//     }
// }

void compute_and_write_dct_scalespace(FILE* fp,
                                      _myfloat* in, int w_in, int h_in,
                                      struct sift_parameters* p,
                                      int flag_interp,
                                      int flag_type)
{
    /* Reading architecture parameters */
    _myfloat sigma_in = p->sigma_in;
    _myfloat sigma_min = p->sigma_min;
    _myfloat delta_min = p->delta_min;
    _myfloat k = pow(2, 1.0/(_myfloat)(p->dog_nspo));
    int nOct = number_of_octaves(w_in, h_in, p);
    int nspo = p->n_spo;

    /* Compute scale-space architecture */
    int* ws  = xmalloc(nOct*sizeof(int));
    int* hs  = xmalloc(nOct*sizeof(int));
    int* nScas = xmalloc(nOct*sizeof(int));
    _myfloat* deltas  = xmalloc(nOct*sizeof(_myfloat));
    _myfloat** sigmas = xmalloc(nOct*sizeof(_myfloat*));
    assert(delta_min <=1);
    deltas[0] = delta_min;
    hs[0] = (int)(h_in/delta_min);
    ws[0] = (int)(w_in/delta_min);
    for(int o = 1; o<nOct; o++){
        ws[o] = ws[o-1]/2;     /*integer division*/
        hs[o] = hs[o-1]/2;
    }
    for(int o = 0; o<nOct; o++){
        nScas[o] = nspo+2;
        deltas[o] = delta_min * pow(2,o); // deltas[o-1]*2.0;
        sigmas[o] = xmalloc(nScas[o]*sizeof(_myfloat));
        for(int s = 0; s < nScas[o]; s++){ /* nSca images + 3 auxiliary images*/
            sigmas[o][s] = sigma_min*pow(2.0, o + (_myfloat)s/(_myfloat)nspo);
        }
    }

    /* Save architecture */
    fwrite(&nOct, sizeof(int), 1, fp);
    for(int o = 0; o < nOct; o++){
        fwrite(&deltas[o], sizeof(_myfloat), 1, fp);
        fwrite(&ws[o], sizeof(int), 1, fp);
        fwrite(&hs[o], sizeof(int), 1, fp);
        fwrite(&nScas[o], sizeof(int), 1, fp);
        fwrite(sigmas[o], sizeof(_myfloat), nScas[o], fp);
    }

    /* Compute interpolated image */
    int w_seed, h_seed;
    _myfloat* seed = compute_interpolated(in, w_in, h_in, p, &w_seed, &h_seed, flag_interp);
    _myfloat* tmp = xmalloc(w_seed*h_seed*sizeof(*tmp)); // same size than the seed


    /* Compute each image of the scalespace */
    for(int o = 0; o < nOct; o++){
        int h = hs[o];
        int w = ws[o];
        int ns = nScas[o];
        int subfactor = pow(2, o);
        for(int s = 0; s < ns; s++){

            // Gaussian convolution + subsampling
            _myfloat sigma = sigmas[o][s] ;
            _myfloat* imM = xmalloc(w*h*sizeof(*imM));
            _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta_min;
            gaussian_blur_with_flag(seed, tmp, w_seed, h_seed, sigM, 1);
            subsample_by_intfactor(tmp, imM, w_seed, h_seed, subfactor);

            if (flag_type == 0){ // Gaussian scale-space
                // Gaussian convolution + subsampling
                _myfloat sigma = sigmas[o][s] ;
                _myfloat* imM = xmalloc(w*h*sizeof(*imM));
                _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta_min;
                gaussian_blur_with_flag(seed, tmp, w_seed, h_seed, sigM, 1);
                subsample_by_intfactor(tmp, imM, w_seed, h_seed, subfactor);
                // SAVE
                fwrite(imM, sizeof(_myfloat), h*w, fp);
            }
            if (flag_type == 1){ // DOG operator
                // image lower scale (in DoG)
                //  Gaussian convolution + subsampling
                _myfloat sigma = sigmas[o][s] ;
                _myfloat* imM = xmalloc(w*h*sizeof(*imM));
                _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta_min;
                gaussian_blur_with_flag(seed, tmp, w_seed, h_seed, sigM, 1);
                subsample_by_intfactor(tmp, imM, w_seed, h_seed, subfactor);
                // image larger scale (in DoG)
                _myfloat* imP = xmalloc(w*h*sizeof(*imP));
                _myfloat sigP = sqrt( k * k * sigma*sigma - sigma_in*sigma_in)/delta_min;
                gaussian_blur_with_flag(seed, tmp, w_seed, h_seed, sigP, 1);
                subsample_by_intfactor(tmp, imP, w_seed, h_seed, subfactor);
                // difference
                for(int i = 0; i < w*h; i++){
                    imM[i] = imP[i] - imM[i];
                }
                // SAVE
                fwrite(imM, sizeof(_myfloat), h*w, fp);
            }
            if (flag_type == 2){
                // Gaussian convolution + subsampling
                _myfloat sigma = sigmas[o][s] ;
                _myfloat* imM = xmalloc(w*h*sizeof(*imM));
                _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta_min;
                gaussian_blur_with_flag(seed, tmp, w_seed, h_seed, sigM, 1);
                subsample_by_intfactor(tmp, imM, w_seed, h_seed, subfactor);
                // Gradient norm + normalization
                _myfloat* im_x = (_myfloat*)xmalloc(h*w*sizeof(_myfloat));
                _myfloat* im_y = (_myfloat*)xmalloc(h*w*sizeof(_myfloat));
                sift_compute_gradient(imM, im_x, im_y, w, h);
                for(int i = 0; i < w*h; i++){
                    imM[i] = sigma*sqrt(im_x[i]*im_x[i] + im_y[i]*im_y[i]);
                }
                // SAVE
                fwrite(imM, sizeof(_myfloat), h*w, fp);
            }
            if (flag_type == 3){
                // Laplacian + normalization
                _myfloat* lapl = (_myfloat*)xmalloc(h*w*sizeof(_myfloat));
                compute_laplacian(imM, lapl, w, h);
                for(int i = 0; i < w*h; i++){
                    imM[i] = sigma*sigma*lapl[i];
                }
                // SAVE
                fwrite(imM, sizeof(_myfloat), h*w, fp);
            
            }
        }
    }
    free(seed);
}



void compute_and_write_dct_scalespace_gradient_norm(FILE* fp,
                                      _myfloat* in, int w_in, int h_in,
                                      struct sift_parameters* p,
                                      int flag_interp,
                                      int flag_type)
{
    /* Reading architecture parameters */
    _myfloat sigma_in = p->sigma_in;
    _myfloat sigma_min = p->sigma_min;
    _myfloat delta_min = p->delta_min;
    _myfloat k = pow(2, 1.0/(_myfloat)(p->dog_nspo));
    int nOct = number_of_octaves(w_in, h_in, p);
    int nspo = p->n_spo;

    /* Compute scale-space architecture */
    int* ws  = xmalloc(nOct*sizeof(int));
    int* hs  = xmalloc(nOct*sizeof(int));
    int* nScas = xmalloc(nOct*sizeof(int));
    _myfloat* deltas  = xmalloc(nOct*sizeof(_myfloat));
    _myfloat** sigmas = xmalloc(nOct*sizeof(_myfloat*));
    assert(delta_min <=1);
    deltas[0] = delta_min;
    hs[0] = (int)(h_in/delta_min);
    ws[0] = (int)(w_in/delta_min);
    for(int o = 1; o<nOct; o++){
        ws[o] = ws[o-1]/2;     /*integer division*/
        hs[o] = hs[o-1]/2;
    }
    for(int o = 0; o<nOct; o++){
        nScas[o] = nspo+2;
        deltas[o] = delta_min * pow(2,o); // deltas[o-1]*2.0;
        sigmas[o] = xmalloc(nScas[o]*sizeof(_myfloat));
        for(int s = 0; s < nScas[o]; s++){ /* nSca images + 3 auxiliary images*/
            sigmas[o][s] = sigma_min*pow(2.0, o + (_myfloat)s/(_myfloat)nspo);
        }
    }

    /* Save architecture */
    fwrite(&nOct, sizeof(int), 1, fp);
    for(int o = 0; o < nOct; o++){
        fwrite(&deltas[o], sizeof(_myfloat), 1, fp);
        fwrite(&ws[o], sizeof(int), 1, fp);
        fwrite(&hs[o], sizeof(int), 1, fp);
        fwrite(&nScas[o], sizeof(int), 1, fp);
        fwrite(sigmas[o], sizeof(_myfloat), nScas[o], fp);
    }

    /* Compute interpolated image */
    int w_seed, h_seed;
    _myfloat* seed = compute_interpolated(in, w_in, h_in, p, &w_seed, &h_seed, flag_interp);
    _myfloat* tmp = xmalloc(w_seed*h_seed*sizeof(*tmp)); // same size than the seed


    /* Compute each image of the scalespace */
    for(int o = 0; o < nOct; o++){
        int h = hs[o];
        int w = ws[o];
        int ns = nScas[o];
        int subfactor = pow(2, o);
        for(int s = 0; s < ns; s++){
            // image at the correct scale
            _myfloat sigma = sigmas[o][s] ;
            _myfloat* imM = xmalloc(w*h*sizeof(*imM));
            _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta_min;
            gaussian_blur_with_flag(seed, tmp, w_seed, h_seed, sigM, 1);
            subsample_by_intfactor(tmp, imM, w_seed, h_seed, subfactor);
            // Gradient norm normalized
            _myfloat* im_x = (_myfloat*)xmalloc(h*w*sizeof(_myfloat));
            _myfloat* im_y = (_myfloat*)xmalloc(h*w*sizeof(_myfloat));
            sift_compute_gradient(imM, im_x, im_y, w, h);
            for(int ii = 0; ii < w*h; ii++){
                im_x[ii] = sigmas[o][s]*sqrt(im_x[ii]*im_x[ii] + im_y[ii]*im_y[ii]);
            }
            fwrite(im_x, sizeof(_myfloat), h*w, fp);
        }
    }
    free(seed);
}






// MEMO TODO this must be suboptimal, 
// TODO TO BE SPLIT IN HALF - to avoid reading the scalespace architecture twice
void crop_from_scalespace_file(_myfloat* cube, FILE* fp, int o, int s, int i, int j, int ri, int rj, int rs, _myfloat ext)
{
    // READ architecture
    int nOct;
    fseek(fp, 0, SEEK_SET); // safety
    fread(&nOct, sizeof(int), 1, fp);
    int* ws  = xmalloc(nOct*sizeof(int));
    int* hs  = xmalloc(nOct*sizeof(int));
    int* nScas  = xmalloc(nOct*sizeof(int));
    _myfloat* deltas  = xmalloc(nOct*sizeof(_myfloat));
    _myfloat** sigmas = xmalloc(nOct*sizeof(_myfloat*));
    for(int i = 0; i < nOct; i++){
        fread(&deltas[i], sizeof(_myfloat), 1, fp);
        fread(&ws[i], sizeof(int), 1, fp);
        fread(&hs[i], sizeof(int), 1, fp);
        fread(&nScas[i], sizeof(int), 1, fp);
        sigmas[i] = xmalloc(nScas[i]*sizeof(_myfloat));
        fread(sigmas[i], sizeof(_myfloat), nScas[i], fp);
        //debug("Architecture Oct=%i (delta, w, h, nScas) = (%f,%i, %i, %i) \n", i, deltas[i], ws[i], hs[i], nScas[i]);
    }

    // stream offset to the first image stack the  of the stream
    int offset;
    fgetpos(fp, &offset);
    for(int i = 0; i < o; i++){
        offset += nScas[i]*ws[i]*hs[i]*sizeof(_myfloat);
    }

    // Extract portion of scalespace
    int h = hs[o];
    int w = ws[o];
    int nSca = nScas[o];
    int n = 0;
    for(int ks = s-rs; ks <= s+rs; ks++){
        for(int ki = i-ri; ki <= i+ri; ki++){
            for(int kj = j-rj; kj <= j+rj; kj++){
                if ( (ks>=0) && (ki>=0) && (kj>=0) && (ks<nSca) && (ki<h) && (kj<w) ){
                    //int ppp;
                    int pos = (ks*w*h + ki*w + kj)*sizeof(_myfloat);
                    fseek(fp, offset+pos, SEEK_SET);
                    fread(&cube[n], sizeof(_myfloat), 1, fp);
                    //fgetpos(fp, &ppp); debug("fread, current position in stream - %i", ppp);
                }else{
                    cube[n] = ext; // possibly NAN
                }
                //fprintf(stderr, "%i %f\n", n, cube[n]);
                n++;
            }
        }
    }

    // Deallocating memory
    xfree(ws);
    xfree(hs);
    xfree(deltas);
    xfree(nScas);
    for(int i = 0; i < nOct; i++){
        // xfree(sigmas[i]);
    }
    xfree(sigmas);
}





