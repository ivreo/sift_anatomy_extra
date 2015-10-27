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

void oversample_with_flag_BIS(const _myfloat* in, int win, int hin, _myfloat* out, int wout, int hout, _myfloat delta_out, int flag_interp);

static void scalespace_compute_dct(struct sift_scalespace* ss,
                            _myfloat* image,
                            int im_w,
                            int im_h,
                            _myfloat sigma_in,
                            int flag_interp)
{

    /* the characteristic of the sead */
    _myfloat delta_min = ss->octaves[0]->delta;
    _myfloat sigma_min = ss->octaves[0]->sigmas[0];
    int w_min = ss->octaves[0]->w;
    int h_min = ss->octaves[0]->h;

    /* checking scale-space definition consistance*/
    assert(w_min == (int)(im_w/delta_min));
    assert(h_min == (int)(im_h/delta_min));

    int nOct = ss->nOct;  /* # of octaves in the scalespace */
    int nSca, w, h;      /* current octave dimensions */
    _myfloat delta;        /* current octave  inter-sample distance */
    _myfloat sig_prev,sig_next,sigma_extra;

    /* for dereferencing */
    struct octa* octave;
    struct octa* octave_prev;
    int w_prev, h_prev;
    _myfloat* im_prev;
    _myfloat* im_next;

    for(int o = 0; o < nOct; o++){

        octave = ss->octaves[o];
        nSca = octave->nSca;
        w = octave->w;
        h = octave->h;
        delta  = octave->delta;  /* intersample distance */

        /** first image in the stack */
        if(o==0){ /* from input image */

            assert(sigma_min>=sigma_in);
            sigma_extra = sqrt(sigma_min*sigma_min - sigma_in*sigma_in)/delta_min;
            if(delta_min < 1){
                _myfloat* tmp = (_myfloat*)malloc(w* h * sizeof(_myfloat));
                oversample_with_flag_BIS(image, im_w, im_h, tmp, w, h, delta, flag_interp);
                add_gaussian_blur_dct(tmp,octave->imStack,w,h,sigma_extra);
                free(tmp);
            }else{ /* ie delta_min = 1, w_min = in_w... */
                add_gaussian_blur_dct(image,octave->imStack,w,h,sigma_extra);
            }
        }
        else{ /* from previous octave */
            octave_prev = ss->octaves[o-1];
            w_prev  = octave_prev->w;
            h_prev = octave_prev->h;
            sift_subsample_by2(&octave_prev->imStack[ (nSca-3) *w_prev*h_prev], octave->imStack, w_prev, h_prev);  /* HiH */
        }

        /** The rest of the image stack*/
        for(int s = 1; s < nSca; s++){ /*add blur to previous image in the stack*/
            im_prev = &octave->imStack[(s-1)*w*h];
            im_next = &octave->imStack[s*w*h];

            sig_prev = octave->sigmas[s-1];
            sig_next = octave->sigmas[s];
            sigma_extra = sqrt(sig_next*sig_next- sig_prev*sig_prev)/delta;
            add_gaussian_blur_dct(im_prev,im_next,w,h,sigma_extra);
        }
    }
}







static void scalespace_compute_nosemigroup(struct sift_scalespace* ss,
                                           _myfloat* image,
                                           int im_w,
                                           int im_h,
                                           _myfloat sigma_in,
                                           int flag_dct,
                                           int flag_interp)
{


    /* the characteristic of the sead */
    _myfloat delta_min = ss->octaves[0]->delta;
    int w_min = ss->octaves[0]->w;
    int h_min = ss->octaves[0]->h;

    /* checking scale-space definition consistance*/
    assert(w_min == (int)(im_w/delta_min));
    assert(h_min == (int)(im_h/delta_min));

    int nOct = ss->nOct;  /* # of octaves in the scalespace */
    _myfloat sigma_abs, sigma_extra;

    // Compute the seed (oversample if necessary)
    int h_seed, w_seed;
    _myfloat delta_seed;
    _myfloat* tmp;
    if(w_min > im_w){
        h_seed = h_min;
        w_seed = w_min;
        delta_seed = delta_min;
        tmp = malloc(w_seed * h_seed * sizeof(_myfloat));
        oversample_with_flag_BIS(image, im_w, im_h, tmp, w_seed, h_seed, delta_seed, flag_interp);
    }else{
        h_seed = im_h;
        w_seed = im_w;
        delta_seed = 1.0;
        tmp = malloc(w_seed * h_seed * sizeof(_myfloat));
        for(int i = 0; i < w_seed * h_seed; i++){
            tmp[i] = image[i];
        }
    }

    for(int o = 0; o < nOct; o++){
        struct octa* octave = ss->octaves[o];
        int ns = octave->nSca;
        int w  = octave->w;
        int h = octave->h;
        for(int s = 0; s < ns; s++){ /*add blur to previous image in the stack*/
            _myfloat *tmp2 = malloc(w_seed* h_seed * sizeof(_myfloat));
            sigma_abs = ss->octaves[o]->sigmas[s];
            sigma_extra = sqrt( sigma_abs*sigma_abs - sigma_in*sigma_in)/delta_seed;
            if (flag_dct == 1){
                add_gaussian_blur_dct( tmp, tmp2, w_seed, h_seed, sigma_extra);
            }else{
                sift_add_gaussian_blur( tmp, tmp2, w_seed, h_seed, sigma_extra);
            }
            int subfactor = pow(2,o);
            subsample_by_intfactor( tmp2, &octave->imStack[s*w*h], w_seed, h_seed, subfactor);
            free(tmp2);
        }
    }
    free(tmp);
}



static void scalespace_compute_LoG(struct sift_scalespace *ss,
                            struct sift_scalespace *nl,
                            int flag_dct)
{


    for(int o = 0; o < nl->nOct; o++){

        struct octa* nl_oct = nl->octaves[o];
        struct octa* ss_oct  = ss->octaves[o];
        int ns = nl_oct->nSca;
        int w = nl_oct->w;
        int h = nl_oct->h;

        for(int s = 0; s < ns; s++){
            /* deferencing */
            _myfloat sigma = ss_oct->sigmas[s];
            _myfloat * lapl = &nl_oct->imStack[s*w*h];
            _myfloat * im   = &ss_oct->imStack[s*w*h];
            if (flag_dct == 1){
                laplacian_with_dct(im, lapl, w, h);
            }
            else{
                compute_laplacian(im, lapl, w, h);
                //compute_laplacian_scheme(im, lapl, w, h, 0.5);
            }
            for(int i = 0; i < h*w; i++){
                lapl[i] *= sigma*sigma;
            }
        }
    }
}

static void scalespace_compute_with_or_without_semigroup(struct sift_scalespace* s,
                                                  _myfloat* x,
                                                  int w,
                                                  int h,
                                                  _myfloat sigma_in,
                                                  int flag_semigroup,
                                                  int flag_dct,
                                                  int flag_interp)
{
    if (flag_semigroup == 1){
        if (flag_dct == 0){
            scalespace_compute(s, x, w, h, sigma_in);
           // scalespace_compute_interpolation(s, x, w, h, sigma_in, flag_interp); // DOESN'T EXIST - TODO
        }
        else{
            scalespace_compute_dct(s, x, w, h, sigma_in, flag_interp);
        }
    }else{
        scalespace_compute_nosemigroup(s, x, w, h, sigma_in, flag_dct, flag_interp);
    }
}




// To avoid descriptor computation
void copy_all_keys(struct sift_keypoints* kA, struct sift_keypoints* kB)
{
    for(int i = 0; i < kA->size; i++){
        struct keypoint* onekey = sift_malloc_keypoint_from_model_and_copy(kA->list[i]);
        sift_add_keypoint_to_list(onekey, kB);
    }
}






_myfloat* scalespace_fiber(int* n, const struct sift_scalespace* ss, _myfloat* xs,  _myfloat* ys, int np, struct sift_parameters* p)
{
    int noct = ss->nOct;
   // int nspo = ss->octaves[0]->nSca - 3; // -2 si DoG;
    int nspo = p->n_spo; // -2 si DoG;
    _myfloat* fiber = malloc(nspo*noct*np*sizeof(_myfloat));
    *n = noct*nspo;
    for(int o = 0; o < noct; o++){
        int w = ss->octaves[o]->w;
        int h = ss->octaves[o]->h;
        _myfloat delta = ss->octaves[o]->delta;
        for (int p = 0;  p < np; p++){
            _myfloat x = xs[p];
            _myfloat y = ys[p];
            int i = (int)(x/delta + 0.5);
            int j = (int)(y/delta + 0.5);
            //for(int s = 0; s < nspo +1 ; s++){
            for(int s = 1; s < nspo +1 ; s++){
                fiber[p*(nspo*noct)+o*nspo+(s-1)] = ss->octaves[o]->imStack[s*w*h+i*w+j];
                //fiber[p*(nspo*noct)+o*nspo+s] = ss->octaves[o]->imStack[s*w*h+i*w+j];
            }
        }
    }
    return fiber;
}





//// used to define the coordinates used to extract a scalespace fiber
//static _myfloat puiss2(_myfloat x)
//{
//    int p = (int)(log(x)/log(2));
//    _myfloat y = pow(2, p);
//    return(y);
//}




static void gaussian_blur_with_flag(_myfloat* in, _myfloat* out, int w, int h, _myfloat sigma, int flag_dct);




void scalespace_compute_dog_controlled(const struct sift_scalespace *s,
                                   struct sift_scalespace *d,
                                   _myfloat k,
                                   int flag_dct)
{
    for(int o = 0; o < d->nOct; o++){
        const struct octa* s_oct = s->octaves[o];
        struct octa* d_oct = d->octaves[o];
        int ns = d_oct->nSca;
        int w = d_oct->w;
        int h = d_oct->h;
        for(int s = 0; s < ns; s++){

            _myfloat* imM = &s_oct->imStack[s*w*h];
            _myfloat* imP = xmalloc(w*h*sizeof(_myfloat));

            _myfloat sigM = d_oct->sigmas[s];
            _myfloat delta = d_oct->delta;
            _myfloat sigP = sigM*sqrt( k*k-1)/delta;
            gaussian_blur_with_flag(imM, imP, w, h, sigP, flag_dct);

            for(int i = 0; i < w*h; i++){
                d_oct->imStack[s*w*h+i] = imP[i] - imM[i];
            }
            free(imP);
        }
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
        //deltas[o] = deltas[o-1]*2.0;
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






struct sift_keypoints* sift_anatomy_dense(_myfloat* x, int w, int h,
                                          struct sift_parameters* p,
                                          struct sift_scalespace* ss[4],
                                          struct sift_keypoints* kk[6],
                                          int flag_semigroup,
                                          int flag_dct,
                                          int flag_log,
                                          int flag_interp)
{
    struct sift_keypoints* k = sift_malloc_keypoints();
    int n_oct = number_of_octaves(w, h, p);

    // normalizing the threshold
    _myfloat thresh = convert_threshold(p);
    if (flag_log){
        //thresh /= exp( M_LN2/( (_myfloat)(p->n_spo))) - 1;
        thresh /=  pow(2, (_myfloat)(p->dog_nspo)) - 1.;
    }

    /** MEMORY ALLOCATION **/
    /** scale-space structure */
    struct sift_scalespace* s   = sift_malloc_scalespace_lowe_floatnspo(n_oct, p->fnspo, w, h, p->delta_min, p->sigma_min);
    struct sift_scalespace* d   = sift_malloc_scalespace_dog_from_scalespace(s);
    struct sift_scalespace* sx  = sift_malloc_scalespace_from_model(s);
    struct sift_scalespace* sy  = sift_malloc_scalespace_from_model(s);
    /** list-of-keypoints (Already allocated) */
    struct sift_keypoints* kA   = kk[0];  /* 3D (discrete) extrema   */
    struct sift_keypoints* kB   = kk[1];  /* passing the threshold on DoG  */
    struct sift_keypoints* kC   = kk[2];  /* interpolated 3D extrema (continuous) */
    struct sift_keypoints* kD   = kk[3];  /* passing the threshold on DoG  */
    struct sift_keypoints* kE   = kk[4];  /* passing OnEdge filter */
    struct sift_keypoints* kF   = kk[5];  /* image border */

    /** KEYPOINT DETECTION ***************************************************/
    scalespace_compute_with_or_without_semigroup(s, x, w, h, p->sigma_in, flag_semigroup, flag_dct, flag_interp);
    if (flag_log == 1){
        scalespace_compute_LoG(s,d, flag_dct);
    }else{
        _myfloat kfact = pow(2, 1.0/(_myfloat)(p->dog_nspo)); // the equivalent factor in the definition of the DoG operator
        scalespace_compute_dog_controlled(s, d, kfact, flag_dct);
    }
    keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
    keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
    keypoints_interpolate_position_controlled(d, kB, kC, p->itermax, p->ofstMax_X, p->ofstMax_S, p->flag_jumpinscale, p->fnspo, p->sigma_min);
    keypoints_discard_with_low_response(kC, kD, thresh);
    keypoints_compute_edge_response(d,kD);
    keypoints_discard_on_edge( kD, k, (p->C_edge+1)*(p->C_edge+1)/p->C_edge );
    (void)kE;
    (void)kF;

    /** scalespace structures*/
    ss[0] = s;
    ss[1] = d;
    ss[2] = sx;
    ss[3] = sy;

    return k;
}




//%//%//%// struct sift_keypoints* sift_anatomy_dense_XXXXXXXXXXXXX(_myfloat* x, int w, int h,
//%//%//%//                                           struct sift_parameters* p,
//%//%//%//                                           struct sift_scalespace* ss[4],
//%//%//%//                                           struct sift_keypoints* kk[6],
//%//%//%//                                           int flag_semigroup,
//%//%//%//                                           int flag_dct,
//%//%//%//                                           int flag_log,
//%//%//%//                                           int flag_interp)
//%//%//%// {
//%//%//%//     struct sift_keypoints* k = sift_malloc_keypoints();
//%//%//%//     int n_oct = number_of_octaves(w, h, p);
//%//%//%// 
//%//%//%//     // normalizing the threshold
//%//%//%//     _myfloat thresh = convert_threshold(p);
//%//%//%//     if (flag_log){
//%//%//%//         //thresh /= exp( M_LN2/( (_myfloat)(p->n_spo))) - 1;
//%//%//%//         thresh /=  pow(2, (_myfloat)(p->dog_nspo)) - 1.;
//%//%//%//     }
//%//%//%// 
//%//%//%//     /** MEMORY ALLOCATION **/
//%//%//%//     /** scale-space structure */
//%//%//%//     struct sift_scalespace* s   = sift_malloc_scalespace_lowe_floatnspo(n_oct, p->fnspo, w, h, p->delta_min, p->sigma_min);
//%//%//%//     struct sift_scalespace* d   = sift_malloc_scalespace_dog_from_scalespace(s);
//%//%//%//     struct sift_scalespace* sx  = sift_malloc_scalespace_from_model(s);
//%//%//%//     struct sift_scalespace* sy  = sift_malloc_scalespace_from_model(s);
//%//%//%//     /** list-of-keypoints (Already allocated) */
//%//%//%//     struct sift_keypoints* kA   = kk[0];  /* 3D (discrete) extrema   */
//%//%//%//     struct sift_keypoints* kB   = kk[1];  /* passing the threshold on DoG  */
//%//%//%//     struct sift_keypoints* kC   = kk[2];  /* interpolated 3D extrema (continuous) */
//%//%//%//     struct sift_keypoints* kD   = kk[3];  /* passing the threshold on DoG  */
//%//%//%//     struct sift_keypoints* kE   = kk[4];  /* passing OnEdge filter */
//%//%//%//     struct sift_keypoints* kF   = kk[5];  /* image border */
//%//%//%// 
//%//%//%// 
//%//%//%// 
//%//%//%//     /** SCALE SPACE COMPUTATION ***************************************************/
//%//%//%//     // Compute the interpolated image (not the seed) - an image with dmin and blur sin (not smin)
//%//%//%//     int w_seed, h_seed;
//%//%//%//     _myfloat* seed = compute_interpolated(x, w, h, p, &w_seed, &h_seed, flag_dct, flag_interp);
//%//%//%// 
//%//%//%// 
//%//%//%// 
//%//%//%// 
//%//%//%// // 
//%//%//%//     scalespace_compute_with_or_without_semigroup(s, x, w, h, p->sigma_in, flag_semigroup, flag_dct, flag_interp);
//%//%//%//     if (flag_log == 1){
//%//%//%//         scalespace_compute_LoG(s,d, flag_dct);
//%//%//%//     }else{
//%//%//%//         _myfloat kfact = pow(2, 1.0/(_myfloat)(p->dog_nspo)); // the equivalent factor in the definition of the DoG operator
//%//%//%//         ///   WARNING - this uses the image at sigma to compute the image at kappa*sigma !!! (gradual is different)
//%//%//%//         scalespace_compute_dog_controlled(s, d, kfact, flag_dct);
//%//%//%//     }
//%//%//%// 
//%//%//%// 
//%//%//%// 
//%//%//%// 
//%//%//%//     /** KEYPOINT DETECTION ***************************************************/
//%//%//%// 
//%//%//%// 
//%//%//%//     keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
//%//%//%//     keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
//%//%//%//     keypoints_interpolate_position_controlled(d, kB, kC, p->itermax, p->ofstMax_X, p->ofstMax_S, p->flag_jumpinscale, p->fnspo, p->sigma_min);
//%//%//%//     keypoints_discard_with_low_response(kC, kD, thresh);
//%//%//%//     keypoints_compute_edge_response(d,kD);
//%//%//%//     keypoints_discard_on_edge( kD, ktemp, (p->C_edge+1)*(p->C_edge+1)/p->C_edge );
//%//%//%//     (void)kE;
//%//%//%//     (void)kF;
//%//%//%// 
//%//%//%//     /** scalespace structures*/
//%//%//%//     ss[0] = s;
//%//%//%//     ss[1] = d;
//%//%//%//     ss[2] = sx;
//%//%//%//     ss[3] = sy;
//%//%//%//     return k;
//%//%//%// }
//%//%//%// 






























static void extract_patch_from_scalespace(_myfloat* list,
                                    struct sift_scalespace* ss,
                                    int i,
                                    int j,
                                    int o,
                                    int s)
{
    // initialisation
    for(int i = 0 ; i < 33*33; i++)
        list[i] = 0;  // sufficient because we print DoG patches


    struct octa* octave = ss->octaves[o];
    int w = octave->w;
    int h = octave->h;
    _myfloat* im = &(ss->octaves[o]->imStack[s*h*w]);
    for (int di = -16; di<=16; di++){
        int ii = (i+di);
        if ((ii>=0)&&(ii<=w-1)){
            for (int dj = -16; dj<=16; dj++){
                int jj = (j+dj);
                if ((jj >= 0)&&(jj <= h-1)){
                    list[(16+di)*33+(16+dj)] = im[ii*w+jj];
                }
            }
        }
    }
}






_myfloat* sift_anatomy_dense_patch(struct sift_keypoints* k,
                                          _myfloat* x, int w, int h,
                                          struct sift_parameters* p,
                                          struct sift_scalespace* ss[4],
                                          struct sift_keypoints* kk[6],
                                          int flag_semigroup,
                                          int flag_dct,
                                          int flag_log,
                                          int flag_interp)
{

    int n_oct = number_of_octaves(w, h, p);



    /** MEMORY ALLOCATION **/
    /** scale-space structure */
    struct sift_scalespace* s   = sift_malloc_scalespace_lowe(n_oct,p->n_spo,w,h,p->delta_min,p->sigma_min);
    struct sift_scalespace* d   = sift_malloc_scalespace_dog_from_scalespace(s);
    /** list-of-keypoints (Already allocated) */
    struct sift_keypoints* kA   = kk[0];  /* 3D (discrete) extrema   */
    struct sift_keypoints* kB   = kk[1];  /* passing the threshold on DoG  */
    struct sift_keypoints* kC   = kk[2];  /* interpolated 3D extrema (continuous) */
    struct sift_keypoints* kD   = kk[3];  /* passing the threshold on DoG  */
    struct sift_keypoints* kE   = kk[4];  /* passing OnEdge filter */
    struct sift_keypoints* kF   = kk[5];  /* image border */

    /** KEYPOINT DETECTION ***************************************************/

    scalespace_compute_with_or_without_semigroup(s, x, w, h, p->sigma_in, flag_semigroup, flag_dct, flag_interp);
    if (flag_log == 1){
        scalespace_compute_LoG(s,d, flag_dct);
    }else{
        scalespace_compute_dog(s,d);
    }

    _myfloat* patch = xmalloc(k->size * 33*33*sizeof(_myfloat));

    for(int id = 0; id < k->size; id++){
        int k_i = k->list[id]->i;
        int k_j = k->list[id]->j;
        int k_o = k->list[id]->o;
        int k_s = k->list[id]->s;
        extract_patch_from_scalespace(&(patch[id*33*33]), d, k_i, k_j, k_o, k_s);
    }

    (void)kA;
    (void)kB;
    (void)kC;
    (void)kD;
    (void)kE;
    (void)kF;

    /** scalespace structures*/
    ss[0] = s;
    ss[1] = d;
   // ss[2] = sx;
   // ss[3] = sy;

   return patch;
}


// Memory allocation for the computation of one octave (with nspo+2 scales)
struct sift_scalespace* sift_malloc_scalespace_dog_gradual_floatnspo_ENTIREOCTAVE(int ObjnOct,
                                                       _myfloat ObjnSca,
                                                       int im_w, int im_h,
                                                       _myfloat delta_min, _myfloat sigma_min,
                                                       int curr_o) // current coordinates
{

    //fprintf(stderr, " In malloc entireoctave  ObjnSca %f\n", ObjnSca );

    int nOct = 1;
    int nSca = (int)(ObjnSca)+2;
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

    for(int s=0;s<nSca;s++){ /* nSca images + 2 auxiliary images*/
        sigmas[0][s] = fsigma_min*pow(2.0,(_myfloat)s / ObjnSca);
        //fprintf(stderr, "------- ------- simulated s %i  scale %33.30f \n", s, sigmas[0][s]);
        //fprintf(stderr, "------- ------- simulated scale %33.30f \n", sigmas[0][s]);
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
        //if (s == 0){
        //fprintf(stderr, "------- ------- simulated  scale %33.30f \n", sigmas[0][s+1]);
        //}
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



//struct sift_scalespace* sift_malloc_scalespace_dog_gradual_floatnspo_auxil(int ObjnOct,
//                                                       _myfloat ObjnSca,
//                                                       int im_w, int im_h,
//                                                       _myfloat delta_min, _myfloat sigma_min,
//                                                       int curr_o, int curr_s) // current coordinates
//{
//    // only compute a octave of three consecutive scales.
//    int nOct = 1;
//    int nSca = 3; //pour le dog
//    // normal
//    int* ws = xmalloc(nOct*sizeof(int));
//    int* hs = xmalloc(nOct*sizeof(int));
//    int* nScas = xmalloc(nOct*sizeof(int));
//    _myfloat* deltas = xmalloc(nOct*sizeof(_myfloat));
//    _myfloat** sigmas = xmalloc(nOct*sizeof(_myfloat*));
//
//    deltas[0] = delta_min * pow(2, curr_o);
//    ws[0] = (int)(im_w / (delta_min * pow(2, curr_o)) );
//    hs[0] = (int)(im_h / (delta_min * pow(2, curr_o)) );
//    _myfloat fsigma_min = sigma_min * pow(2, curr_o);
//    nScas[0] = nSca;
//    sigmas[0] = xmalloc(nScas[0]*sizeof(_myfloat));
//    for(int s = -1; s <= 1; s++){
//        sigmas[0][s+1] = fsigma_min*pow(2.0,(_myfloat)(s+curr_s)/ObjnSca);
//    }
//    struct sift_scalespace* scalespace = sift_malloc_scalespace(nOct, deltas, ws, hs, nScas, sigmas);
//    xfree(deltas);
//    xfree(ws);
//    xfree(hs);
//    xfree(nScas);
//    xfree(sigmas[0]);
//    xfree(sigmas);
//    return scalespace;
//}
//
//
//

























// TODO const _myfloat* in, but lib_fourier complains
static void gaussian_blur_with_flag(_myfloat* in, _myfloat* out, int w, int h, _myfloat sigma, int flag_dct)
{
    if (flag_dct == 1){
        add_gaussian_blur_dct(in, out, w, h, sigma);
    }
    else{
        sift_add_gaussian_blur(in, out, w, h, sigma);
    }
}


/**
 *    @ oversampling with various interpolation methods
 *    flag_interp values:
 *      0 : bilinear
 *      1 : DCT
 *      3,5,7,9,11 : Bsplines
 *
 */
void oversample_with_flag_BIS(const _myfloat* in, int win, int hin,
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






_myfloat* compute_interpolated(const _myfloat* in, int w, int h, struct sift_parameters* p, int* w_seed, int* h_seed, int flag_dct, int flag_interp)
{
    _myfloat delta = p->delta_min;
    int wout = (int)(w/delta);
    int hout = (int)(h/delta);
    _myfloat* out = malloc(wout*hout*sizeof(*out));
    // oversampling
    oversample_with_flag_BIS(in, w, h, out, wout, hout, delta, flag_interp);
    *w_seed = wout;
    *h_seed = hout;
    return(out);
}



// Computes a one octave scale-space
//   - images are computed from the interpolated image

void compute_controlled_nosemi_dog_from_interpolated(_myfloat* in, int w_in, int h_in, _myfloat delta_in, _myfloat sigma_in,
                                                 struct sift_scalespace *d, _myfloat k, int flag_dct, int subfactor)
{
    _myfloat* tmp = xmalloc(w_in*h_in*sizeof(*tmp));
    struct octa* d_oct = d->octaves[0];
    int ns = d_oct->nSca;
    int w = d_oct->w;
    int h = d_oct->h;
    //
    for(int s = 0; s < ns; s++){

        // image lower scale (in DoG)
        _myfloat sigma = d_oct->sigmas[s] ;
        _myfloat* imM = xmalloc(w*h*sizeof(*imM));
        _myfloat sigM = sqrt(sigma*sigma - sigma_in*sigma_in)/delta_in;
        gaussian_blur_with_flag(in, tmp, w_in, h_in, sigM, flag_dct);
        subsample_by_intfactor(tmp, imM, w_in, h_in, subfactor);
        //fprintf(stderr, "computing  s=%i / sigma %f / sigma supp %f / k %f \n", s, sigma, sigM, k);

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



void update_controlled_nosemi_dog_from_interpolated(_myfloat* in, int w_in, int h_in, _myfloat delta_in, _myfloat sigma_in,
                                                 struct sift_scalespace *d, _myfloat k, int flag_dct, int subfactor)
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









void keypoints_discard_with_scale_above(struct sift_keypoints *keysIn,
                                        struct sift_keypoints *keysAccept,
                                        _myfloat max_scale)
{
    for( int k = 0; k < keysIn->size; k++){
        struct keypoint* key = keysIn->list[k];
        bool isAccepted = ( key->sigma <  max_scale );
        if (isAccepted == true){
            struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy, keysAccept);
        }
    }
}




_myfloat convert_threshold_DoG(const struct sift_parameters* p) // A copy of convert_threshold (but uses dog_nspo instead of n_spo)
{
    // converting the threshold to make it consistent
    //_myfloat k_nspo =  exp( M_LN2/( (_myfloat)(p->dog_nspo)));
    //_myfloat k_3 =  exp( M_LN2/( (_myfloat)3));
    _myfloat k_nspo = pow(2, (_myfloat)(p->dog_nspo));
    _myfloat k_3 =  pow(2, 3.0);

    _myfloat thresh = (k_nspo - 1) / (k_3 - 1) * p->C_DoG ;
    return thresh;
}



/*
 *
 *  First compute the seed image and then compute the scalespace 
 *
 *
 */
struct sift_keypoints* sift_anatomy_gradual(_myfloat* x, int w, int h,
                                            struct sift_parameters* p,
                                            struct sift_scalespace* ss[4],
                                            struct sift_keypoints* kk[6],
                                            int flag_semigroup,
                                            int flag_dct,
                                            int flag_log,
                                            int flag_interp)
{
    struct sift_keypoints* k = sift_malloc_keypoints();

    int n_oct = number_of_octaves(w, h, p);

    _myfloat thresh;
    if (flag_log){
        thresh = convert_threshold(p);
        thresh /=  pow(2, (_myfloat)(p->dog_nspo)) - 1.;
    }
    else{ // DoG then
        thresh = convert_threshold_DoG(p);
    }

    // WARNING - HERE we call seed the interpolated image WARNING its blur is sigma_in
    int w_seed, h_seed;
    _myfloat* seed = compute_interpolated(x, w, h, p, &w_seed, &h_seed, flag_dct, flag_interp);



   // for(int o = 0; o < n_oct; o++){
    for(int o = 0; o < n_oct; o++){

        struct sift_scalespace* d = sift_malloc_scalespace_dog_gradual_floatnspo(n_oct, p->fnspo ,w,h, p->delta_min,p->sigma_min, o, 1);
        int subfactor = pow(2, o);
        int nscales = ceil(p->fnspo);
        for(int is = 0; is < nscales; is++){ // current scale where the keypoints are detected

            /** list-of-keypoints (Already allocated) */
            struct sift_keypoints* kA = sift_malloc_keypoints();
            struct sift_keypoints* kB = sift_malloc_keypoints();
            struct sift_keypoints* kC = sift_malloc_keypoints();
            struct sift_keypoints* kD = sift_malloc_keypoints();
            struct sift_keypoints* kE = sift_malloc_keypoints();
            //
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
        
                   // if (s == 0){
                   //     fprintf(stderr, "------- ------- simulated scale %33.30f \n", d->octaves[0]->sigmas[s]);
                   // }

                }
                // the missing image
                update_controlled_nosemi_dog_from_interpolated(seed, w_seed, h_seed, p->delta_min, p->sigma_in,
                                                                   d, factDoG, flag_dct, subfactor);
            }

       //     //// Note Mauricio : plus rapide d'ecrire en tiff 
       //     //// // DoG scalespace
       //     char name[FILENAME_MAX];
       //     sprintf(name, "grad_oct%04i_scale%04i",  o, is+1); // warning - decalage indice
       //     print_sift_scalespace_rgb(d, name);
       //     //print_sift_scalespace_gray(d, name);

            if (p->itermax == 0){  // All discrete extrema
                keypoints_find_3d_discrete_extrema_epsilon(d, ktmp, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
            }else{
                keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
                keypoints_interpolate_position_controlled(d, kA, ktmp, p->itermax, p->ofstMax_X, p->ofstMax_S, 0, p->fnspo, pow(2,o)*p->sigma_min);
                debug(" current scale = %i  / Number of detections: discrete=%i  interpolated=%i", is, kA->size, ktmp->size);

               // // BEGIN COMMENT Following is the the normal flow of keypoint.
               // keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
               // keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
               // keypoints_interpolate_position_controlled(d, kB, kC, p->itermax, p->ofstMax_X, p->ofstMax_S, 0, p->fnspo, pow(2,o)*p->sigma_min);
               // keypoints_discard_with_low_response(kC, kD, thresh);
               // keypoints_compute_edge_response(d,kD);
               // keypoints_discard_on_edge(kD, ktmp, (p->C_edge+1)*(p->C_edge+1)/p->C_edge);
               // // END COMMENT.
               // //


            }

            // 20150612
            // update the index of octave and scale (by default all keypoints are in scale 1 (gradual implementation)).
            for(int ik = 0; ik < ktmp->size; ik++){
                struct keypoint* key = ktmp->list[ik];
                key->o = o;
                key->s = is+1;
                struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
                sift_add_keypoint_to_list(copy, k);
                //fprintf(stderr, "DEBUG o = %i / is = %i  /////  k last element  o %i  s %i  \n ", o, is, k->list[k->size -1]->o, k->list[k->size -1]->s  );
            }

            sift_free_keypoints(kA);
            sift_free_keypoints(kB);
            sift_free_keypoints(kC);
            sift_free_keypoints(kD);
            sift_free_keypoints(kE);
            //
            sift_free_keypoints(ktmp);
        }
        sift_free_scalespace(d);
    }
    free(seed);
    return k;
}



// 
//
static void keypoints_only_in_given_octave(struct sift_keypoints *keysIn,
                                           struct sift_keypoints *keysAccept,
                                           int o)
{
    for( int k = 0; k < keysIn->size; k++){
        struct keypoint *key = keysIn->list[k];
        bool isAccepted = (key->o == o);
        if (isAccepted == true){
            struct keypoint *copy = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy, keysAccept);
        }
    }
}


//   based on gradual - the entire octaves are computed sequentially
void sift_anatomy_gradual_check_extrema(_myfloat* x, int w, int h,
                                        struct sift_keypoints* keysIn,
                                        struct sift_keypoints* keysTested, // keypoints from keysIn that can be tested.
                                        struct sift_keypoints* keysExtrema,
                                        struct sift_keypoints* keysNotExtrema,
                                        struct sift_parameters* p,
                                        int flag_semigroup,
                                        int flag_dct,
                                        int flag_interp)
{
    int n_oct = number_of_octaves(w, h, p);

    _myfloat thresh  = convert_threshold_DoG(p);

    // Computation of the interpolated image
    int w_seed, h_seed;
    _myfloat* seed = compute_interpolated(x, w, h, p, &w_seed, &h_seed, flag_dct, flag_interp);

    for(int o = 0; o < n_oct; o++){

        // considering only the keypoints that are in the current octave
        struct sift_keypoints* keysOct = sift_malloc_keypoints();
        keypoints_only_in_given_octave(keysIn, keysOct, o);
        struct sift_keypoints* keysOctExtrema = sift_malloc_keypoints();
        struct sift_keypoints* keysOctNotExtrema = sift_malloc_keypoints();

        struct sift_scalespace* d = sift_malloc_scalespace_dog_gradual_floatnspo_ENTIREOCTAVE(n_oct, p->fnspo ,w,h, p->delta_min,p->sigma_min, o);
        int subfactor = pow(2, o);


        // Compute single octave scale-space
        _myfloat factDoG = pow(2, 1.0/(_myfloat)(p->dog_nspo));
        compute_controlled_nosemi_dog_from_interpolated(seed, w_seed, h_seed, p->delta_min, p->sigma_in,
                                                                d, factDoG, flag_dct, subfactor);

        // Check the extrema
        //keypoints_check_3d_discrete_extrema_epsilon_sphere(d, keysOct, keysOctExtrema, keysOctNotExtrema, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->discrete_sphere_r, p->discrete_sphere_dr);
        //keypoints_confirm_extremum_present_in_ball(d, keysOct, keysOctExtrema, keysOctNotExtrema, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->ball_alpha, p->ball_beta);
        keypoints_confirm_extremum_present_in_ball(d, keysOct, keysOctExtrema, keysOctNotExtrema, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->ball_r1, p->ball_r2, p->ball_r3);

        // Update the keypoint octave index (by default all keypoints are in octave 0 (gradual implementation)).
        for(int ik = 0; ik < keysOctExtrema->size; ik++){
            struct keypoint* key = keysOctExtrema->list[ik];
            // copy the keypoint to the list of confirmed extrema
            struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy, keysExtrema);
            // copy the keypoint to the list of tested keypoints
            struct keypoint* copy2 = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy2, keysTested);
        }
        for(int ik = 0; ik < keysOctNotExtrema->size; ik++){
            struct keypoint* key = keysOctNotExtrema->list[ik];
            // copy the keypoint to the list of keys not confirmed as extrema
            struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy, keysNotExtrema);
            // copy the keypoint to the list of tested keypoints
            struct keypoint* copy2 = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy2, keysTested);
        }


        // 
        sift_free_keypoints(keysOctNotExtrema);
        sift_free_keypoints(keysOctExtrema);
        sift_free_keypoints(keysOct);
        sift_free_scalespace(d);
    }
    free(seed);
}



//   
struct sift_keypoints* sift_anatomy_gradual_ENTIREOCTAVE(_myfloat* x, int w, int h,
                                            struct sift_parameters* p,
                                            struct sift_scalespace* ss[4],
                                            struct sift_keypoints* kk[6],
                                            int flag_semigroup,
                                            int flag_dct,
                                            int flag_log,
                                            int flag_interp)
{
    struct sift_keypoints* k = sift_malloc_keypoints();

    int n_oct = number_of_octaves(w, h, p);

    _myfloat thresh;
    if (flag_log){
        thresh = convert_threshold(p);
        thresh /=  pow(2, (_myfloat)(p->dog_nspo)) - 1.;
    }
    else{ // DoG then
        thresh = convert_threshold_DoG(p);
    }

    // Computation of the interpolated image
    int w_seed, h_seed;
    _myfloat* seed = compute_interpolated(x, w, h, p, &w_seed, &h_seed, flag_dct, flag_interp);


    for(int o = 0; o < n_oct; o++){
        
        //struct sift_scalespace* d = sift_malloc_scalespace_dog_gradual_floatnspo(n_oct, p->fnspo ,w,h, p->delta_min,p->sigma_min, o, 1);
        struct sift_scalespace* d = sift_malloc_scalespace_dog_gradual_floatnspo_ENTIREOCTAVE(n_oct, p->fnspo ,w,h, p->delta_min,p->sigma_min, o);
        int subfactor = pow(2, o);

        /** list-of-keypoints (Already allocated) */
        struct sift_keypoints* kA = sift_malloc_keypoints();
        struct sift_keypoints* kB = sift_malloc_keypoints();
        struct sift_keypoints* kC = sift_malloc_keypoints();
        struct sift_keypoints* kD = sift_malloc_keypoints();
        struct sift_keypoints* kE = sift_malloc_keypoints();
        //
        struct sift_keypoints* ktmp = sift_malloc_keypoints(); // to update the octave and scale indices
        


        struct sift_keypoints* kExtrema = sift_malloc_keypoints();    // to update the octave and scale indices
        struct sift_keypoints* kNotExtrema = sift_malloc_keypoints(); // to update the octave and scale indices


        /** KEYPOINT DETECTION ***************************************************/
        _myfloat factDoG = pow(2, 1.0/(_myfloat)(p->dog_nspo)); // the equivalent factor in the definition of the DoG operator
        compute_controlled_nosemi_dog_from_interpolated(seed, w_seed, h_seed, p->delta_min, p->sigma_in,
                                                                d, factDoG, flag_dct, subfactor);

 //       //// Note Mauricio : plus rapide d'ecrire en tiff 
 //       //// // DoG scalespace
 //       char name[FILENAME_MAX];
 //       sprintf(name, "grad_oct%04i_",  o);
 //       print_sift_scalespace_rgb(d, name);
 //       //print_sift_scalespace_gray(d, name);

        if (p->itermax == 0){  // All discrete extrema
            //keypoints_find_3d_discrete_extrema_epsilon(d, ktmp, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
            //keypoints_find_3d_discrete_extrema_epsilon_largervolume(d, ktmp, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->discrete_extrema_dR);

            keypoints_find_3d_discrete_extrema_epsilon_sphere(d, ktmp, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->discrete_sphere_r, p->discrete_sphere_dr);


 // //       // TEST - the check routine
 // //           
 // //           // keypoints_find_3d_discrete_extrema_epsilon_sphere(d, ktmp, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->discrete_sphere_r, p->discrete_sphere_dr);
 // //           keypoints_find_3d_discrete_extrema_epsilon(d, ktmp, p->n_ori, p->n_hist, p->n_bins, p->epsilon);

 // //            keypoints_check_3d_discrete_extrema_epsilon_sphere(d, ktmp, kExtrema, kNotExtrema, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->discrete_sphere_r, p->discrete_sphere_dr);

 // //            fprintf(stderr, "test check sphere  nextrema:%i  checktrue:%i  checkfalse%i \n", ktmp->size, kExtrema->size, kNotExtrema->size);    




        }else{
            //fprintf(stderr, "before extraction \n");
            //keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
            //keypoints_find_3d_discrete_extrema_epsilon_largervolume(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->discrete_extrema_dR);
            keypoints_find_3d_discrete_extrema_epsilon_sphere(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon, p->discrete_sphere_r, p->discrete_sphere_dr);
            //fprintf(stderr, "after extraction \n");
            keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
            keypoints_interpolate_position_controlled(d, kB, kC, p->itermax, p->ofstMax_X, p->ofstMax_S, 0, p->fnspo, pow(2,o)*p->sigma_min);
            keypoints_discard_with_low_response(kC, kD, thresh);
            keypoints_compute_edge_response(d,kD);
            keypoints_discard_on_edge(kD, ktmp, (p->C_edge+1)*(p->C_edge+1)/p->C_edge);
        }

        // 20150612
        // update the keypoint octave index (by default all keypoints are in octave 0 (gradual implementation)).
        for(int ik = 0; ik < ktmp->size; ik++){
            struct keypoint* key = ktmp->list[ik];
            key->o = o;
            struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy, k);
        }


        sift_free_keypoints(kA);
        sift_free_keypoints(kB);
        sift_free_keypoints(kC);
        sift_free_keypoints(kD);
        sift_free_keypoints(kE);
        //
        sift_free_keypoints(ktmp);
        sift_free_scalespace(d);
        
        
        
        //
        sift_free_keypoints(kExtrema);
        sift_free_keypoints(kNotExtrema);
    }
    free(seed);
    return k;
}










struct sift_keypoints* sift_anatomy_gradual_auxil(_myfloat* x, int w, int h,
                                            struct sift_parameters* p,
                                            struct sift_scalespace* ss[4],
                                            struct sift_keypoints* kk[6],
                                            int flag_semigroup,
                                            int flag_dct,
                                            int flag_log,
                                            int flag_interp)
{
    struct sift_keypoints* k = sift_malloc_keypoints();

    int n_oct = number_of_octaves(w, h, p);

    _myfloat thresh;
    if (flag_log){
        thresh = convert_threshold(p);
        thresh /= pow(2, (_myfloat)(p->dog_nspo)) - 1.;
    }
    else{ // DoG then
        thresh = convert_threshold_DoG(p);
    }

    // WARNING - HERE we call seed the interpolated image WARNING its blur is sigma_in
    int w_seed, h_seed;
    _myfloat* seed = compute_interpolated(x, w, h, p, &w_seed, &h_seed, flag_dct, flag_interp);
    

   // for(int o = 0; o < n_oct; o++){
    for(int o = 0; o < n_oct; o++){

        //struct sift_scalespace* d = sift_malloc_scalespace_dog_gradual_floatnspo( n_oct, p->fnspo ,w,h, p->delta_min,p->sigma_min, o, 1);
        struct sift_scalespace* d = sift_malloc_scalespace_dog_gradual_floatnspo( n_oct, p->fnspo ,w,h, p->delta_min,p->sigma_min, o, 0);
        // CHECK that smin* 2^(-1/nspo) > sin  (for extra auxiliary scale)
        assert(d->octaves[0]->sigmas[0] > p->sigma_in );

        int subfactor = pow(2, o);
        int nscales = ceil(p->fnspo);
        // auxiliary scale
        for(int is = -1; is < nscales; is++){ // TEMP AUXILIARY - one extra scale

            /** list-of-keypoints (Already allocated) */
            struct sift_keypoints* kA = sift_malloc_keypoints();
            struct sift_keypoints* kB = sift_malloc_keypoints();
            struct sift_keypoints* kC = sift_malloc_keypoints();
            struct sift_keypoints* kD = sift_malloc_keypoints();
            struct sift_keypoints* kE = sift_malloc_keypoints();

            /** KEYPOINT DETECTION ***************************************************/
            _myfloat factDoG = pow(2, 1.0/(_myfloat)(p->dog_nspo)); // the equivalent factor in the definition of the DoG operator
            if( is == -1){
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
    //        //// Note Mauricio : plus rapide d'ecrire en tiff 
    //        //// // DoG scalespace
    //        char name[FILENAME_MAX];
    //        sprintf(name, "grad_auxil_oct%04i_scale%04i", o, 1+is); // TEMP decalage d'indice.
    //        //print_sift_scalespace_rgb(d, name);
    //        print_sift_scalespace_gray(d, name);

            if (p->itermax == 0){  // All discrete extrema
                keypoints_find_3d_discrete_extrema_epsilon(d, k, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
            }else{
                keypoints_find_3d_discrete_extrema_epsilon(d, kA, p->n_ori, p->n_hist, p->n_bins, p->epsilon);
                keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
                keypoints_interpolate_position_controlled(d, kB, kC, p->itermax, p->ofstMax_X, p->ofstMax_S, 0, p->fnspo, pow(2,o)*p->sigma_min);
                keypoints_discard_with_low_response(kC, kD, thresh);
                keypoints_compute_edge_response(d,kD);
                keypoints_discard_on_edge(kD, k, (p->C_edge+1)*(p->C_edge+1)/p->C_edge);
            }

            sift_free_keypoints(kA);
            sift_free_keypoints(kB);
            sift_free_keypoints(kC);
            sift_free_keypoints(kD);
            sift_free_keypoints(kE);
        }
        sift_free_scalespace(d);
    }
    free(seed);
    return k;
}
















//   /* @brief compute the seed image (first image in the scalespace (sigma_min))
//    *
//    *
//    */
//   _myfloat* compute_seed(const _myfloat* in, int w, int h, struct sift_parameters* p, int* w_seed, int* h_seed, int flag_dct, int flag_interp)
//   {
//       _myfloat delta = p->delta_min;
//       _myfloat sigin = p->sigma_in;
//       _myfloat sigmin = p->sigma_min;
//       int wout = (int)(w/delta);
//       int hout = (int)(h/delta);
//       _myfloat* tmp = malloc(wout*hout*sizeof(*tmp));
//       _myfloat* out = malloc(wout*hout*sizeof(*out));
//       // oversampling
//       oversample_with_flag_BIS(in, w, h, tmp, wout, hout, delta, flag_interp);
//       // extra blur
//       _myfloat sigma = sqrt(sigmin*sigmin - sigin*sigin)/delta;
//       gaussian_blur_with_flag(tmp, out, wout, hout, sigma, flag_dct);
//       // output
//       *w_seed = wout;
//       *h_seed = hout;
//       free(tmp);
//       return(out);
//   
//   }
//
//
//
//
//   /*
//    *
//    *   The definition of the DoG operator is controlled independently from the scalespace grid sampling.
//    *   @param k
//    *
//    */
//   void scalespace_compute_dog_from_interpolated_controlled(_myfloat* seed, int wseed, int hseed, _myfloat deltaseed, _myfloat sigmaseed,
//                                                    struct sift_scalespace *d, _myfloat k, int flag_dct, int subfactor)
//   {
//       _myfloat* tmp = malloc(wseed*hseed*sizeof(*tmp));
//       //for(int o = 0; o < d->nOct; o++){ // TEMP
//       for(int o = 0; o < 1; o++){
//           fprintf(stderr, "    kkkkkkkk  %i sur %i \n", o, d->nOct);
//           struct octa* d_oct = d->octaves[o];
//           int ns = d_oct->nSca;
//           int w = d_oct->w;
//           int h = d_oct->h;
//           for(int s = 0; s < ns; s++){
//               _myfloat sigma = d_oct->sigmas[s];
//               _myfloat* imM = malloc(w*h*sizeof(*imM));
//               _myfloat sigM = sqrt(sigma*sigma - sigmaseed*sigmaseed)/deltaseed;
//               gaussian_blur_with_flag(seed, tmp, wseed, hseed, sigM, flag_dct);
//               subsample_by_intfactor(tmp, imM, wseed, hseed, subfactor);
//               // note : using semi-group imP computed from imM
//               _myfloat delta = d_oct->delta;
//               _myfloat sigP = sigma*sqrt(k*k-1.0)/delta;
//               _myfloat* imP = malloc(w*h*sizeof(*imP));
//               gaussian_blur_with_flag(imM, imP, w, h, sigP, flag_dct);
//               for(int i = 0; i < w*h; i++){
//                   d_oct->imStack[s*w*h+i] = imP[i] - imM[i];
//               }
//               free(imP);
//               free(imM);
//           }
//       }
//       free(tmp);
//   }
//   
//   void scalespace_update_dog_from_interpolated_controlled(_myfloat* seed, int wseed, int hseed, _myfloat deltaseed, _myfloat sigmaseed,
//                                                    struct sift_scalespace *d, _myfloat k, int flag_dct, int subfactor)
//   {
//       _myfloat* tmp = malloc(wseed*hseed*sizeof(*tmp));
//       for(int o = 0; o < 1; o++){
//           struct octa* d_oct = d->octaves[o];
//           int ns = d_oct->nSca;
//           int w = d_oct->w;
//           int h = d_oct->h;
//           // NOTE :  in general ns = 3 (gradual computation).
//           // COPY previously computed images
//           for(int i = 0; i < (ns-1)*w*h; i++){
//               d_oct->imStack[i] = d_oct->imStack[i+w*h];
//           }
//           // ADD the new DoG image on top
//           //    - compute the image at scale sigma_new 
//           _myfloat new_sigma = d_oct->sigmas[2];  // supposing it has been correctly updated before
//           _myfloat* imM = malloc(w*h*sizeof(*imM));
//           _myfloat sigM = sqrt(new_sigma*new_sigma - sigmaseed*sigmaseed)/deltaseed;
//           gaussian_blur_with_flag(seed, tmp, wseed, hseed, sigM, flag_dct);
//           //subsample_by_intfactor(tmp, imM, wseed, hseed, (int)( (_myfloat)hseed/(_myfloat)h + 0.5));
//           subsample_by_intfactor(tmp, imM, wseed, hseed, subfactor);
//           fprintf(stderr, "comp subsampling factor   %i  sigM = %.40f  \n", subfactor, sigM);
//           //    - compute the image at scale k*sigma_new 
//           //      (note, this using semi-group)
//           _myfloat delta = d_oct->delta;
//           _myfloat sigP = new_sigma*sqrt(k*k-1.0)/delta;
//           _myfloat* imP = malloc(w*h*sizeof(*imP));
//           gaussian_blur_with_flag(imM, imP, w, h, sigP, flag_dct);
//           //    - compute the difference
//           for(int i = 0; i < w*h; i++){
//               d_oct->imStack[(ns-1)*w*h+i] = imP[i] - imM[i];
//           }
//           free(imP);
//           free(imM);
//       }
//       free(tmp);
//   }



















//// OBSOLETE
//struct sift_scalespace* sift_malloc_scalespace_dog_gradual(int ObjnOct,
//                                                       int ObjnSca,
//                                                       int im_w, int im_h,
//                                                       _myfloat delta_min, _myfloat sigma_min,
//                                                       int curr_o, int curr_s) // current coordinates
//{
//
//    // only compute a octave of three consecutive scales.
//    int nOct = 1;
//    int nSca = 3; //pour le dog
//    // normal
//    int* ws = xmalloc(nOct*sizeof(int));
//    int* hs = xmalloc(nOct*sizeof(int));
//    int* nScas = xmalloc(nOct*sizeof(int));
//    _myfloat* deltas = xmalloc(nOct*sizeof(_myfloat));
//    _myfloat** sigmas = xmalloc(nOct*sizeof(_myfloat*));
//
//    deltas[0] = delta_min;
//    ws[0] = (int)(im_w/delta_min);
//    hs[0] = (int)(im_h/delta_min);
//    for(int o = 1; o <= curr_o; o++){
//        ws[0] /= 2;
//        hs[0] /= 2;
//        deltas[0] *= 2.0;
//    }
//    _myfloat fsigma_min = sigma_min * deltas[0] / delta_min;
//    nScas[0] = nSca;
//    sigmas[0] = xmalloc(nScas[0]*sizeof(_myfloat));
//    for(int s = -1; s <= 1; s++){
//        sigmas[0][s+1] = fsigma_min*pow(2.0,(_myfloat)(s+curr_s)/(_myfloat)ObjnSca);
//    }
//    struct sift_scalespace* scalespace = sift_malloc_scalespace(nOct, deltas, ws, hs, nScas, sigmas);
//    xfree(deltas);
//    xfree(ws);
//    xfree(hs);
//    xfree(nScas);
//    xfree(sigmas[0]);
//    xfree(sigmas);
//    return scalespace;
//
//}


















//static void scalespace_compute_with_or_without_semigroup(struct sift_scalespace* s,
//                                                  _myfloat* x,
//                                                  int w,
//                                                  int h,
//                                                  _myfloat sigma_in,
//                                                  int flag_semigroup,
//                                                  int flag_dct)
//{
//    if (flag_semigroup == 1){
//        if (flag_dct == 0){
//            scalespace_compute(s, x, w, h, sigma_in);
//        }
//        else{
//            scalespace_compute_dct(s, x, w, h, sigma_in);
//        }
//    }else{
//        scalespace_compute_nosemigroup(s, x, w, h, sigma_in, flag_dct);
//    }
//}


