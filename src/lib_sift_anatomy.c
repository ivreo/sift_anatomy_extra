/*
IPOL SIFT
Copyright (C) 2014, Ives Rey-Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20140911 (September 11th, 2014)

This C ANSI source code is related to the IPOL publication

    [1] "Anatomy of the SIFT Method."
        I. Rey Otero  and  M. Delbracio
        Image Processing Online, 2013.
        http://www.ipol.im/pub/algo/rd_anatomy_sift/

An IPOL demo is available at
        http://www.ipol.im/pub/demo/rd_anatomy_sift/



== Patent Warning and License =================================================


This program is free software: you can use, modify and/or
redistribute it under the terms of the simplified BSD
License. You should have received a copy of this license along
this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.



The SIFT method is patented

    [2] "Method and apparatus for identifying scale invariant features
      in an image."
        David G. Lowe
        Patent number: 6711293
        Filing date: Mar 6, 2000
        Issue date: Mar 23, 2004
        Application number: 09/519,89

 These source codes are made available for the exclusive aim of serving as
 scientific tool to verify the soundness and completeness of the algorithm
 description. Compilation, execution and redistribution of this file may
 violate patents rights in certain countries. The situation being different
 for every country and changing over time, it is your responsibility to
 determine which patent rights restrictions apply to you before you compile,
 use, modify, or redistribute this file. A patent lawyer is qualified to make
 this determination. If and only if they don't conflict with any patent terms,
 you can benefit from the following license terms attached to this file.

*/
/**
 * @file lib_sift_anatomy.c
 * @brief SIFT anatomy interface.
 *
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */


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

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )
#define ABS(x) ((x)<0?-(x):(x))
//#define EPSILON 0
       //#define EPSILON 0.00001 // for quad
       //#define EPSILON 0.001 // for quad
#define EPSILON FLT_EPSILON
//#define EPSILON DBL_EPSILON



/** @brief Compute the Gaussian scalespace for SIFT
 *
 *
 *  in @param input image  of im_w X im_h
 *     @param sigma_in : assumed level of blur in the image
 *
 *
 *  The construction parameters are already stored in scalespace
 *
 */
void scalespace_compute(struct sift_scalespace* ss,
                        const _myfloat* image,
                        int im_w,
                        int im_h,
                        _myfloat sigma_in)
//
{
    // seed image
    _myfloat delta_min = ss->octaves[0]->delta;
    _myfloat sigma_min = ss->octaves[0]->sigmas[0];
    int w_min = ss->octaves[0]->w;
    int h_min = ss->octaves[0]->h;
    // check that the scalespace is correctly defined
    assert(w_min == (int)(im_w/delta_min));
    assert(h_min == (int)(im_h/delta_min));

    _myfloat sig_prev,sig_next,sigma_extra;

    /* for dereferencing */
    struct octa* octave;
    struct octa* octave_prev;
    int w_prev, h_prev;
    _myfloat* im_prev;
    _myfloat* im_next;

    int nOct = ss->nOct;
    for(int o = 0; o < nOct; o++){

        octave = ss->octaves[o];
        int nSca = octave->nSca;
        int w = octave->w;
        int h = octave->h;
        _myfloat delta  = octave->delta;  /* intersample distance */

        /** first image in the stack */
        if(o==0){ /* from input image */
            assert(sigma_min>=sigma_in);
            sigma_extra = sqrt(sigma_min*sigma_min - sigma_in*sigma_in)/delta_min;
            if(delta_min < 1){
                _myfloat* imtmp = xmalloc(w* h * sizeof(_myfloat));
                sift_oversample_bilin(image,im_w,im_h,imtmp,w,h,delta_min);
                sift_add_gaussian_blur(imtmp,octave->imStack,w,h,sigma_extra);
                xfree(imtmp);
            }else{ /* ie delta_min = 1, w_min = in_w... */
                sift_add_gaussian_blur(image,octave->imStack,w,h,sigma_extra);
            }
        }
        else{ /* from previous octave */
            octave_prev = ss->octaves[o-1];
            w_prev  = octave_prev->w;
            h_prev = octave_prev->h;
            // WARNING nSca also includes the three auximliary pictures.
            int nspo = (nSca-3);  // scales per octace
            sift_subsample_by2(&octave_prev->imStack[nspo*w_prev*h_prev], octave->imStack, w_prev, h_prev);
        }

        /** The rest of the image stack*/
        for(int s = 1; s < nSca; s++){ /*add blur to previous image in the stack*/
            im_prev = &octave->imStack[(s-1)*w*h];
            im_next = &octave->imStack[s*w*h];
            sig_prev = octave->sigmas[s-1];
            sig_next = octave->sigmas[s];
            sigma_extra = sqrt(sig_next*sig_next- sig_prev*sig_prev)/delta;
            sift_add_gaussian_blur(im_prev,im_next,w,h,sigma_extra);
        }
    }
}




/** @brief difference of Gaussians
 *
 */
//static void scalespace_compute_dog(const struct sift_scalespace *s, // TEMP (for lib_dense_anatomy.c)
void scalespace_compute_dog(const struct sift_scalespace *s,
                                   struct sift_scalespace *d)
{
    for(int o = 0; o < d->nOct; o++){
        const struct octa* s_oct = s->octaves[o];
        struct octa* d_oct = d->octaves[o];
        int ns = d_oct->nSca;
        int w = d_oct->w;
        int h = d_oct->h;
        for(int s = 0; s < ns; s++){
            _myfloat* diff = &d_oct->imStack[s*w*h];
            _myfloat* im_P = &s_oct->imStack[(s+1)*w*h];
            _myfloat* im_M = &s_oct->imStack[s*w*h];
            for(int p=0;p<w*h;p++)
                diff[p]=im_P[p]-im_M[p];
        }
    }
}


/** @brief Compute the 2d gradient of each image in the scale-space
 *
 *  @in scalespace
 *  @out dx_scalespace,
 *       scalespace structure storing the gradient x-component of each image
 *
 *  @out dy_scalespace,
 *       scalespace structure storing the gradient y-component of each image
 *
 *
 * The gradients of the auxiliary images are not computed.
 *
 */
void scalespace_compute_gradient(const struct sift_scalespace* scalespace,
                                 struct sift_scalespace* sx,
                                 struct sift_scalespace* sy)
{
    int nOct = scalespace->nOct;
    for(int o = 0; o < nOct; o++){
        int nSca   = scalespace->octaves[o]->nSca; //WARNING this includes the auxiliary images.
        int w  = scalespace->octaves[o]->w;
        int h = scalespace->octaves[o]->h;
        for(int s = 0;s < nSca;s++){
            const _myfloat* im = &scalespace->octaves[o]->imStack[s*w*h];
            _myfloat* dx = &sx->octaves[o]->imStack[s*w*h];
            _myfloat* dy = &sy->octaves[o]->imStack[s*w*h];
            sift_compute_gradient(im, dx, dy, w, h);
        }
    }
}


/** @brief Extract discrete extrema from DoG scale-space.
 *
 * in  @param dog  : Difference of Gaussian (scale-space structure)
 *
 * out @param keys : List of keypoints.
 *
 *
 *
 *    The following parameters are for memory allocation only.
 *
 * @param n_bins : Number of bins in the histogram
 *                 used to attribute principal orientation.
 *
 * @param n_hist
 * @param n_ori  : The SIFT descriptor is an array of
 *                 n_hist X n_hist weighted orientation
 *                 with n_ori bins each.
 *
 */
//STATIC
void keypoints_find_3d_discrete_extrema(struct sift_scalespace* d,
                                               struct sift_keypoints* keys,
                                               int n_ori,
                                               int n_hist,
                                               int n_bins)
{
    for(int o = 0; o < d->nOct; o++){

        int ns = d->octaves[o]->nSca; // dimension of the image stack in the current octave
        int w = d->octaves[o]->w;
        int h = d->octaves[o]->h;
        _myfloat delta = d->octaves[o]->delta; // intersample distance
        _myfloat* imStack = d->octaves[o]->imStack;

        // Precompute index offsets to the 26 neighbors in a 3x3x3 neighborhood.
        int neighbor_offsets[26];
        int n = 0;
        for (int ds = -1; ds <= 1; ds++) {
            for (int di = -1; di <= 1; di++) {
                for (int dj = -1; dj <= 1; dj++) {
                    if (ds != 0 || di != 0 || dj != 0) {
                        neighbor_offsets[n] = (ds * h + di) * w + dj;
                        n++;
                    }
                }
            }
        }

        // EXTRA - the 27 values (the center and its 26 neighbors on the 26 grid
        int neighbor_offsets27[27];
        n = 0;
        for (int ds = -1; ds <= 1; ds++) {
            for (int di = -1; di <= 1; di++) {
                for (int dj = -1; dj <= 1; dj++) {
                    neighbor_offsets27[n] = (ds * h + di) * w + dj;
                    n++;
                }
            }
        }

        // Loop through the samples of the image stack (one octave)
        for(int s = 1; s < ns-1; s++){
            for(int i = 1; i < h-1; i++){
                for(int j = 1; j < w-1; j++){

                    const _myfloat* center = &imStack[s*w*h+i*w+j];
                    const _myfloat center_value = *center;

                    bool is_local_min = true;
                    // An optimizing compiler will unroll this loop.
                    for (int n = 0; n < 26; n++) {
                        if (center[neighbor_offsets[n]] - EPSILON <= center_value) {
                            is_local_min = false;
                            break; // Can stop early if a smaller neighbor was found.
                        }
                    }
                    bool is_local_max = true;
                    // Can skip max check if center point was determined to be a local min.
                    if (is_local_min) {
                        is_local_max = false;
                    } else {
                        // An optimizing compiler will unroll this loop.
                        for (int n = 0; n < 26; n++) {
                            if (center[neighbor_offsets[n]] + EPSILON >= center_value) {
                                is_local_max = false;
                                break; // Can stop early if a larger neighbor was found.
                            }
                        }
                    }
                    if(is_local_max || is_local_min) { // if 3d discrete extrema, save a candidate keypoint
                        struct keypoint* key = sift_malloc_keypoint(n_ori, n_hist, n_bins);
                        key->i = i;
                        key->j = j;
                        key->s = s;
                        key->o = o;
                        key->x = delta*i;
                        key->y = delta*j;
                        key->sigma = d->octaves[o]->sigmas[s];
                        key->val = imStack[s*w*h+i*w+j];
                        //j++; // skip the next pixel (it's not an extremum)
                        // EXTRA - the 27 values (the center and its 26 neighbors on the 26 grid
                        for(int in = 0; in < 27; in++){
                            key->neighbors[in] = imStack[s*w*h+i*w+j + neighbor_offsets27[in]];
                        }
                        sift_add_keypoint_to_list(key,keys);
                    }
                }
            }
        }
    }
}





/////////////////    EPSILON    ////////////////////////////
/////////////////    EPSILON    ////////////////////////////
/////////////////    EPSILON    ////////////////////////////
/////////////////    EPSILON    ////////////////////////////
void keypoints_find_3d_discrete_extrema_epsilon(struct sift_scalespace* d,
                                               struct sift_keypoints* keys,
                                               int n_ori,
                                               int n_hist,
                                               int n_bins,
                                               _myfloat epsilon)
{
    for(int o = 0; o < d->nOct; o++){

        int ns = d->octaves[o]->nSca; // dimension of the image stack in the current octave

        int w = d->octaves[o]->w;
        int h = d->octaves[o]->h;
        _myfloat delta = d->octaves[o]->delta; // intersample distance
        _myfloat* imStack = d->octaves[o]->imStack;

        // Precompute index offsets to the 26 neighbors in a 3x3x3 neighborhood.
        int neighbor_offsets[26];
        int n = 0;
        for (int ds = -1; ds <= 1; ds++) {
            for (int di = -1; di <= 1; di++) {
                for (int dj = -1; dj <= 1; dj++) {
                    if (ds != 0 || di != 0 || dj != 0) {
                        neighbor_offsets[n] = (ds * h + di) * w + dj;
                        n++;
                    }
                }
            }
        }
        
        // EXTRA - the 27 values (the center and its 26 neighbors on the 26 grid
        int neighbor_offsets27[27];
        n = 0;
        for (int ds = -1; ds <= 1; ds++) {
            for (int di = -1; di <= 1; di++) {
                for (int dj = -1; dj <= 1; dj++) {
                    neighbor_offsets27[n] = (ds * h + di) * w + dj;
                    n++;
                }
            }
        }

        // Loop through the samples of the image stack (one octave)
        for(int s = 1; s < ns-1; s++){
            for(int i = 1; i < h-1; i++){
                for(int j = 1; j < w-1; j++){

                    const _myfloat* center = &imStack[s*w*h+i*w+j];
                    const _myfloat center_value = *center;

                    bool is_local_min = true;
                    // An optimizing compiler will unroll this loop.
                    for (int n = 0; n < 26; n++) {
                        if (center[neighbor_offsets[n]] - epsilon <= center_value) {
                            is_local_min = false;
                            break; // Can stop early if a smaller neighbor was found.
                        }
                    }
                    bool is_local_max = true;
                    // Can skip max check if center point was determined to be a local min.
                    if (is_local_min) {
                        is_local_max = false;
                    } else {
                        // An optimizing compiler will unroll this loop.
                        for (int n = 0; n < 26; n++) {
                            if (center[neighbor_offsets[n]] + epsilon >= center_value) {
                                is_local_max = false;
                                break; // Can stop early if a larger neighbor was found.
                            }
                        }
                    }
                    if(is_local_max || is_local_min) { // if 3d discrete extrema, save a candidate keypoint
                        struct keypoint* key = sift_malloc_keypoint(n_ori, n_hist, n_bins);
                        key->i = i;
                        key->j = j;
                        key->s = s;
                        key->o = o;
                        key->x = delta*i;
                        key->y = delta*j;
                        key->sigma = d->octaves[o]->sigmas[s];
                        key->val = imStack[s*w*h+i*w+j];
                        // keep the 27 values (the center and its 26 neighbors on the 26 grid
                        for(int in = 0; in < 27; in++){
                            key->neighbors[in] = imStack[s*w*h+i*w+j + neighbor_offsets27[in]];
                        }
                        sift_add_keypoint_to_list(key,keys);
                    }
                }
            }
        }
    }
}





void keypoints_find_3d_discrete_extrema_epsilon_largervolume(struct sift_scalespace* d,
                                               struct sift_keypoints* keys,
                                               int n_ori,
                                               int n_hist,
                                               int n_bins,
                                               _myfloat epsilon,
                                               int dR)  // The sample is compared to the sample on the surface of a (2 x halfR + 1)^3 volume.
{
    for(int o = 0; o < d->nOct; o++){

        int ns = d->octaves[o]->nSca; // dimension of the image stack in the current octave
        int w = d->octaves[o]->w;
        int h = d->octaves[o]->h;
        _myfloat delta = d->octaves[o]->delta; // intersample distance
        _myfloat* imStack = d->octaves[o]->imStack;

        // TODO normalize

        // Precompute index offsets for the samples on the surface of the cube.
        int nComp = 24*dR*dR +2; // number of comparisons // number of elements on the surface.
        int neighbor_offsets[nComp];
        int n = 0;
        for (int ds = -dR; ds <= dR; ds++) {
            for (int di = -dR; di <= dR; di++) {
                for (int dj = -dR; dj <= dR; dj++) {
                    if (ds == -dR || ds == dR || di == -dR || di == dR ||  dj == -dR || dj == dR) {
                        neighbor_offsets[n] = (ds * h + di) * w + dj;
                        //fprintf(stderr, " ---  n %i / offset %i \n ", n, neighbor_offsets[n] );
                        n++;
                    }
                }
            }
        }


        // Loop through the samples of the image stack (one octave)
        for(int s = dR; s < ns-dR; s++){
            for(int i = dR; i < h-dR; i++){
                for(int j = dR; j < w-dR; j++){

                    //fprintf(stderr, " --- loop through the samples in the image stack %i - %i - %i\n", s, i, j);

                    const _myfloat* center = &imStack[s*w*h+i*w+j];
                    const _myfloat center_value = *center;

                    bool is_local_min = true;
                    // An optimizing compiler will unroll this loop.
                    for (int n = 0; n < nComp; n++) {
                        if (center[neighbor_offsets[n]] - epsilon <= center_value) {
                            is_local_min = false;
                            break; // Can stop early if a smaller neighbor was found.
                        }
                    }
                    bool is_local_max = true;
                    // Can skip max check if center point was determined to be a local min.
                    if (is_local_min) {
                        is_local_max = false;
                    } else {
                        // An optimizing compiler will unroll this loop.
                        for (int n = 0; n < nComp; n++) {
                            if (center[neighbor_offsets[n]] + epsilon >= center_value) {
                                is_local_max = false;
                                break; // Can stop early if a larger neighbor was found.
                            }
                        }
                    }
                    if(is_local_max || is_local_min) { // if 3d discrete extrema, save a candidate keypoint
                        struct keypoint* key = sift_malloc_keypoint(n_ori, n_hist, n_bins);
                        key->i = i;
                        key->j = j;
                        key->s = s;
                        key->o = o;
                        key->x = delta*i;
                        key->y = delta*j;
                        key->sigma = d->octaves[o]->sigmas[s];
                        key->val = imStack[s*w*h+i*w+j];
                        //
                        // j++; // skip the next pixel (it's not an extremum)
                        // EXTRA - the 27 values (the center and its 26 neighbors on the 26 grid
                        // TODO there must more than 27 values here
                        // 
                        //      
                    //    for(int in = 0; in < 27; in++){
                    //        key->neighbors[in] = imStack[s*w*h+i*w+j + neighbor_offsets27[in]];
                    //    }
                        for(int in = 0; in < 26; in++){
                            //key->neighbors[in] = imStack[s*w*h+i*w+j + neighbor_offsets27[in]];
                            key->neighbors[in] = imStack[s*w*h+i*w+j + neighbor_offsets[in]];
                        }

                        sift_add_keypoint_to_list(key,keys);
                    }
                }
            }
        }
    }
}



void keypoints_find_3d_discrete_extrema_epsilon_sphere(struct sift_scalespace* d,
                                               struct sift_keypoints* keys,
                                               int n_ori,
                                               int n_hist,
                                               int n_bins,
                                               _myfloat epsilon,
                                               _myfloat r,
                                               _myfloat dr)  // The sample is compared to the sample on the surface of a (2 x halfR + 1)^3 volume.
{
    for(int o = 0; o < d->nOct; o++){

        int ns = d->octaves[o]->nSca; // dimension of the image stack in the current octave
        int w = d->octaves[o]->w;
        int h = d->octaves[o]->h;
        _myfloat delta = d->octaves[o]->delta; // intersample distance
        _myfloat* imStack = d->octaves[o]->imStack;

        // Precompute index offsets for the samples on the surface of the cube.
        int dR = (int)ceil(r+dr);
        int H = 2*dR + 1;
        int nInVol = H*H*H; // number of samples considered in the volume
        int neighbor_offsets[nInVol];
        int n = 0;
        for (int ds = -dR; ds <= dR; ds++) {
            for (int di = -dR; di <= dR; di++) {
                for (int dj = -dR; dj <= dR; dj++) {
                    // compute the sample distance to the candidate extrema
                    _myfloat dist = sqrt( di*di + dj*dj + ds*ds );
                    if ( (dist >= r) && (dist <= (r+dr) )) {
                    //if (ds == -dR || ds == dR || di == -dR || di == dR ||  dj == -dR || dj == dR) {
                        neighbor_offsets[n] = (ds * h + di) * w + dj;
                        //fprintf(stderr, " ---  n %i / offset %i \n ", n, neighbor_offsets[n] );
                        n++;
                    }
                }
            }
        }
        int nComp = n; // number of comparisons
        fprintf(stderr, "Number of comparisons = %i\n", nComp);


        // Loop through the samples of the image stack (one octave)
        for(int s = dR; s < ns-dR; s++){
            for(int i = dR; i < h-dR; i++){
                for(int j = dR; j < w-dR; j++){

                    //fprintf(stderr, " --- loop through the samples in the image stack %i - %i - %i\n", s, i, j);

                    const _myfloat* center = &imStack[s*w*h+i*w+j];
                    const _myfloat center_value = *center;

                    bool is_local_min = true;
                    // An optimizing compiler will unroll this loop.
                    for (int n = 0; n < nComp; n++) {
                        if (center[neighbor_offsets[n]] - epsilon <= center_value) {
                            is_local_min = false;
                            break; // Can stop early if a smaller neighbor was found.
                        }
                    }
                    bool is_local_max = true;
                    // Can skip max check if center point was determined to be a local min.
                    if (is_local_min) {
                        is_local_max = false;
                    } else {
                        // An optimizing compiler will unroll this loop.
                        for (int n = 0; n < nComp; n++) {
                            if (center[neighbor_offsets[n]] + epsilon >= center_value) {
                                is_local_max = false;
                                break; // Can stop early if a larger neighbor was found.
                            }
                        }
                    }
                    if(is_local_max || is_local_min) { // if 3d discrete extrema, save a candidate keypoint
                        struct keypoint* key = sift_malloc_keypoint(n_ori, n_hist, n_bins);
                        key->i = i;
                        key->j = j;
                        key->s = s;
                        key->o = o;
                        key->x = delta*i;
                        key->y = delta*j;
                        key->sigma = d->octaves[o]->sigmas[s];
                        key->val = imStack[s*w*h+i*w+j];
                        //
                        // j++; // skip the next pixel (it's not an extremum)
                        // EXTRA - the 27 values (the center and its 26 neighbors on the 26 grid
                        // TODO there must more than 27 values here
                        // 
                        //      
                    //    for(int in = 0; in < 27; in++){
                    //        key->neighbors[in] = imStack[s*w*h+i*w+j + neighbor_offsets27[in]];
                    //    }
                        for(int in = 0; in < 26; in++){
                            //key->neighbors[in] = imStack[s*w*h+i*w+j + neighbor_offsets27[in]];
                            key->neighbors[in] = imStack[s*w*h+i*w+j + neighbor_offsets[in]];
                        }

                        sift_add_keypoint_to_list(key,keys);
                    }
                }
            }
        }
    }
}






// We check a set of keypoints by comparing the nearest discrete sample on the
// surface of a sphere defined by parameter r and dr.
void keypoints_check_3d_discrete_extrema_epsilon_sphere(struct sift_scalespace* d,
                                               struct sift_keypoints* keysIn,
                                               struct sift_keypoints* keysExtrema,
                                               struct sift_keypoints* keysNotExtrema,
                                               int n_ori,
                                               int n_hist,
                                               int n_bins,
                                               _myfloat epsilon,
                                               _myfloat r,
                                               _myfloat dr)
{


    int nComp=0; // number of value comparison to confirm that the point is an extrema

    // load keypoints to check
    for( int k = 0; k < keysIn->size; k++){

        struct keypoint* key = keysIn->list[k];
       // int o = key->o;
        int s = key->s;
        int i = key->i;
        int j = key->j;


        // Read the octave dimension -  note As far as 'd' is concerned, we are in the first octave.
        int ns = d->octaves[0]->nSca; // dimension of the image stack in the current octave
        int w = d->octaves[0]->w;
        int h = d->octaves[0]->h;
        _myfloat* imStack = d->octaves[0]->imStack;


        // Precompute index offsets for the samples on the surface of the cube.
        //
        int dR = (int)ceil(r+dr);
        int H = 2*dR + 1;
        int nInVol = H*H*H; // number of samples considered in the volume
        int neighbor_offsets[nInVol];
        int n = 0;
        for (int ds = -dR; ds <= dR; ds++) {
            for (int di = -dR; di <= dR; di++) {
                for (int dj = -dR; dj <= dR; dj++) {
                    // compute the sample distance to the candidate extrema
                    _myfloat dist = sqrt( di*di + dj*dj + ds*ds );
                    if ( (dist >= r) && (dist <= (r+dr) )) {
                        neighbor_offsets[n] = (ds * h + di) * w + dj;
                        n++;
                    }
                }
            }
        }
        nComp = n; // number of comparisons
        




        // Test the extrema by comparing the nearest discrete sample to the
        // samples on the surface of the sphere.
        //
        if ( s>=dR && s<ns-dR && i>=dR && i<h-dR && j>=dR && j<w-dR) {


            const _myfloat* center = &imStack[s*w*h+i*w+j];
            const _myfloat center_value = *center;

            bool is_local_min = true;
            // An optimizing compiler will unroll this loop.
            for (int n = 0; n < nComp; n++) {
                if (center[neighbor_offsets[n]] - epsilon <= center_value) {
                    is_local_min = false;
                    break; // Can stop early if a smaller neighbor was found.
                }
            }
            bool is_local_max = true;
            // Can skip max check if center point was determined to be a local min.
            if (is_local_min) {
                is_local_max = false;
            } else {
                // An optimizing compiler will unroll this loop.
                for (int n = 0; n < nComp; n++) {
                    if (center[neighbor_offsets[n]] + epsilon >= center_value) {
                        is_local_max = false;
                        break; // Can stop early if a larger neighbor was found.
                    }
                }
            }
            if(is_local_max || is_local_min) { // if 3d discrete extrema, save a candidate keypoint
                struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
                sift_add_keypoint_to_list(copy, keysExtrema);
            }else{
                struct keypoint* copyOut = sift_malloc_keypoint_from_model_and_copy(key);
                sift_add_keypoint_to_list(copyOut, keysNotExtrema);
            }
        }
    }
    fprintf(stderr, "number of comparisons to confirm it's an extrema %i \n", nComp);
}






// We check a set of keypoints by the scalespace values in the middle of the
// ball to values that are on the surface of the ball
void keypoints_confirm_extremum_present_in_ball(struct sift_scalespace* d,
                                                struct sift_keypoints* keysIn,
                                                struct sift_keypoints* keysExtrema,
                                                struct sift_keypoints* keysNotExtrema,
                                                int n_ori,
                                                int n_hist,
                                                int n_bins,
                                                _myfloat epsilon,
                                                _myfloat r1,
                                                _myfloat r2,
                                                _myfloat r3)
{

    // Read the octave dimension
    // THIS IS A TRICK: As far as 'd' is concerned, we are in the first octave.
    // THIS IS A TRICK: All keypoints are in the current octave, i.e., the only octave actually in memory, i.e., the octave with index 0.
    int ns = d->octaves[0]->nSca; // dimension of the image stack in the current octave
    int w  = d->octaves[0]->w;
    int h  = d->octaves[0]->h;
    _myfloat* imStack = d->octaves[0]->imStack;
    // THIS IS A TRICK: Now we can precompute the neighbors' positions for all keypoints

//    assert( (alpha < beta ) && (beta < 1) );
//    _myfloat r3 =  (ns-2) / 10.0; // corresponding to the height of the volume considered by Lowe
//                                  // (in logscale)
//    _myfloat r1 = alpha*r3;
//    _myfloat r2 = beta*r3;

    FILE* fin = fopen("map_inside.txt", "w");
    FILE* fout = fopen("map_outside.txt", "w");
    FILE* fall = fopen("map_all.txt", "w");

    // Precompute index offsets for the samples on the surface of the cube.
    int dR = (int)ceil(r3);
    int H = 2*dR + 1;
    fprintf(stderr,  "CONFIRMATION: ball dimension -  r1 %f  - r2 %f - r3 %f    H  %i  nspo %i \n", r1, r2, r3,  H, ns);
    int n_cube = H*H*H; // number of samples in the large cuve the ball is in
    //
    int* neighbor_in = xmalloc(n_cube*sizeof(int));
    int* neighbor_out = xmalloc(n_cube*sizeof(int));
    int* neighbor_all = xmalloc(n_cube*sizeof(int));
    //
    int n_in = 0;
    int n_out = 0;
    int n_all = 0;
    for (int ds = -dR; ds <= dR; ds++) {
        for (int di = -dR; di <= dR; di++) {
            for (int dj = -dR; dj <= dR; dj++) {
                // compute the sample distance to the candidate extrema
                _myfloat dist2 = di*di + dj*dj + ds*ds;
                // list of neighbors that are in the middle of the  ball
                if (dist2 < r1*r1) {
                    neighbor_in[n_in] = (ds * h + di) * w + dj;
                    n_in++;
                    // map
                    fprintf(fin, "%i %i %i \n", ds, di, dj);
                }
                // list of neighbors that are on the surface of the ball
                if ((dist2 >= r2*r2) && (dist2 < r3*r3)) {
                    neighbor_out[n_out] = (ds * h + di) * w + dj;
                    n_out++;
                    // map
                    fprintf(fout, "%i %i %i \n", ds, di, dj);
                }
                neighbor_all[n_all] = (ds * h + di) * w + dj;
                n_all++;
                fprintf(fall, "%i %i %i \n", ds, di, dj);
            }
        }
    }
    fprintf(stderr, "DEBUG after setting up the neighbors positions\n");
    assert( (n_in > 0) && (n_out > 0));

    fclose(fout);
    fclose(fin);
    fclose(fall);

    // Arrays with the scalespace values surrounding a keypoint.
    _myfloat values_in[n_in];     // - inside ball
    _myfloat values_out[n_out];   // - on the ball surface
    _myfloat values_all[n_all];   // - the entire cube

    // Same thing (scalespace values in the neighborhood) for all keypoints.
    _myfloat* all_values_in  = xmalloc(n_in*  (keysIn->size) *sizeof(_myfloat));
    _myfloat* all_values_out = xmalloc(n_out* (keysIn->size) *sizeof(_myfloat));
    _myfloat* all_values_all = xmalloc(n_all* (keysIn->size) *sizeof(_myfloat));

    int n_test = -1;  // number of keypoints actually tested (not discarded because of being too close to the border).
    for( int k = 0; k < keysIn->size; k++){

        struct keypoint* key = keysIn->list[k];
        //int o = key->o;
        int s = key->s;
        int i = key->i;
        int j = key->j;
        // DEBUG
        //fprintf(stderr, "DEBUGin_confirm_extrema: Candidate keypoint (s,i,j) = (%i,%i,%i) \n", s,i,j);

        // Confirm that an extremum is present on this point
        if ( s>=dR && s<ns-dR && i>=dR && i<h-dR && j>=dR && j<w-dR){  // some points are inevitably discarded on the octave borders.

            n_test +=1;

            // Load values, compute mins, maxs, means and stds.
            const _myfloat* center = &imStack[s*w*h+i*w+j];
            for (int n = 0; n < n_in; n++){
                values_in[n] = center[neighbor_in[n]];
            }
            _myfloat max_in =  array_max(values_in, n_in);
            _myfloat min_in =  array_min(values_in, n_in);
            _myfloat mean_in, std_in;
            array_mean_and_std(values_in, n_in, &mean_in, &std_in);

            for (int n = 0; n < n_out; n++){
                values_out[n] = center[neighbor_out[n]];
            }
            _myfloat max_out =  array_max(values_out, n_out);
            _myfloat min_out =  array_min(values_out, n_out);
            _myfloat mean_out, std_out;
            array_mean_and_std(values_out, n_out, &mean_out, &std_out);

            // for visualisation only - not used in criteria
            for (int n = 0; n < n_all; n++){
                values_all[n] = center[neighbor_all[n]];
            }

            // criterion
            //bool is_local_max = (min_in > max_out);
            //bool is_local_min = (max_in < min_out);

            // denoised criterion (this is not exactly Mauricio's since the values on the surface are not smoothed locally).
            //bool is_local_max = (mean_in > max_out);
            //bool is_local_min = (mean_in < min_out);

            // JM criterion
            bool is_local_max = (max_in > max_out);
            bool is_local_min = (min_in < min_out);

            if(is_local_max || is_local_min) { // saving the keypoint in the appropriate list of keypoint
                struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
                sift_add_keypoint_to_list(copy, keysExtrema);
            }else{
                struct keypoint* copyOut = sift_malloc_keypoint_from_model_and_copy(key);
                sift_add_keypoint_to_list(copyOut, keysNotExtrema);
            }

            // Record the surrounding for all keypoints.
            // - in the middle.
            for(int n = 0; n < n_in; n++){
                all_values_in[n_test*n_in+n] = values_in[n];
            }
            // - on the surface of the ball.
            for(int n = 0; n < n_out; n++){
                all_values_out[n_test*n_out+n] = values_out[n];
            }
            // - in the entire surrounding cube.
            for(int n = 0; n < n_all; n++){
                all_values_all[n_test*n_all+n] = values_all[n];
            }
        }
    }
    fprintf(stderr, "Number of samples used in the criteria:\n   in the middle:%i\n   on the surface: %i \n", n_in, n_out);


    // BEGIN SAVE DATA
    char name[1024*1024];
    FILE* f;
    // save in binary file
    sprintf(name, "all_values_inside.bin");
    f = fopen(name, "wb");
    fwrite(all_values_in, sizeof(_myfloat), n_test*n_in , f);
    fclose(f);

    sprintf(name, "all_values_outside.bin");
    f = fopen(name, "wb");
    fwrite(all_values_out, sizeof(_myfloat), n_test*n_out , f);
    fclose(f);

    // save in ascii file
    // - in
    sprintf(name, "all_values_inside.txt");
    f = fopen(name, "w");
    for (int k = 0; k< n_test; k++){
        for (int n = 0; n < n_in; n++){
            fprintf(f, "%32.30f ", all_values_in[k*n_in+n]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    // - out
    sprintf(name, "all_values_outside.txt");
    f = fopen(name, "w");
    for (int k = 0; k< n_test; k++){
        for (int n = 0; n < n_out; n++){
            fprintf(f, "%32.30f ", all_values_out[k*n_out+n]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    // - all
    sprintf(name, "all_values_all.txt");
    f = fopen(name, "w");
    for (int k = 0; k< n_test; k++){
        for (int n = 0; n < n_all; n++){
            fprintf(f, "%32.30f ", all_values_all[k*n_all+n]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    // END SAVE DATA

    free(all_values_in);
    free(all_values_out);
    free(all_values_all);

    free(neighbor_in);
    free(neighbor_out);
    free(neighbor_all);
}































/** @brief Filter a list of keypoints (copy the keypoints that satisfy a criteria)
 *
 * in  @param keysIn : list of keys to be filtered.
 * 
 * out @param keysAccept   list of keys with DoG absolute value beyond threshold.
 * 
 *  @param thresh    :constant threshold over the entire scale-space.
 * 
 */
//static void keypoints_discard_with_low_response(struct sift_keypoints *keysIn, // TEMP (lib_dense_anatomy)
void keypoints_discard_with_low_response(struct sift_keypoints *keysIn, // TEMP (lib_dense_anatomy)
                                                struct sift_keypoints *keysAccept,
                                                _myfloat thresh)
{
    for( int k = 0; k < keysIn->size; k++){
        struct keypoint* key = keysIn->list[k];
        bool isAccepted = ( ABS(key->val) >  thresh);
        if (isAccepted == true){
            struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy, keysAccept);
        }
    }
}




/** @brief Execute one interpolation (to refine extrema position)
 * 
 * in  @param imStck      : a stack of images in the DoG space
 *                    w X h X nscales samples.
 *                   
 * in  @param i , j ,  s  : 3d discrete extrema
 * 
 * out @param di , dj , ds : offset in each spatial direction
 * out @param dVal         , the extremum value is : imStck[i,j,s]+dVal
 * 
 * 
 * Computes the 3D Hessian of the DoG space in one point
 * 
 * 
 */
static void inverse_3D_Taylor_second_order_expansion(_myfloat *stack,
                                              int w, int h, int ns,
                                              int i, int j, int s,
                                              _myfloat *di, _myfloat *dj, _myfloat *ds, _myfloat *val)
{
    _myfloat hXX,hXY,hXS,hYY,hYS,hSS;
    _myfloat det,aa,ab,ac,bb,bc,cc;
    _myfloat gX,gY,gS;
    _myfloat ofstX, ofstY, ofstS, ofstVal;
    
    /** Compute the 3d Hessian at pixel (i,j,s)  Finite difference scheme  *****/
    hXX = stack[s*w*h+(i-1)*w+j] + stack[s*w*h+(i+1)*w+j] - 2*stack[s*w*h+i*w+j];
    hYY = stack[s*w*h+i*w+(j+1)] + stack[s*w*h+i*w+(j-1)] - 2*stack[s*w*h+i*w+j];
    hSS = stack[(s+1)*w*h+i*w+j] + stack[(s-1)*w*h+i*w+j] - 2*stack[s*w*h+i*w+j];
    hXY = 0.25*(  (stack[s*w*h+(i+1)*w+(j+1)] - stack[s*w*h+(i+1)*w+(j-1)])
                - (stack[s*w*h+(i-1)*w+(j+1)] - stack[s*w*h+(i-1)*w+(j-1)]) );
    hXS = 0.25*(  (stack[(s+1)*w*h+(i+1)*w+j] - stack[(s+1)*w*h+(i-1)*w+j])
                - (stack[(s-1)*w*h+(i+1)*w+j] - stack[(s-1)*w*h+(i-1)*w+j]) );
    hYS = 0.25*(  (stack[(s+1)*w*h+i*w+(j+1)] - stack[(s+1)*w*h+i*w+(j-1)])
                - (stack[(s-1)*w*h+i*w+(j+1)] - stack[(s-1)*w*h+i*w+(j-1)]) );
    
    /** Compute the 3d gradient at pixel (i,j,s) */
    gX = 0.5*( stack[s*w*h+(i+1)*w+j] - stack[s*w*h+(i-1)*w+j] );
    gY = 0.5*( stack[s*w*h+i*w+(j+1)] - stack[s*w*h+i*w+(j-1)] );
    gS = 0.5*( stack[(s+1)*w*h+i*w+j] - stack[(s-1)*w*h+i*w+j] );
    
    /** Inverse the Hessian - Fitting a quadratic function */
    det = hXX*hYY*hSS - hXX*hYS*hYS - hXY*hXY*hSS + 2*hXY*hXS*hYS - hXS*hXS*hYY ; 
    aa = (hYY*hSS - hYS*hYS)/det;
    ab = (hXS*hYS - hXY*hSS)/det;
    ac = (hXY*hYS - hXS*hYY)/det;
    bb = (hXX*hSS - hXS*hXS)/det;
    bc = (hXY*hXS - hXX*hYS)/det;
    cc = (hXX*hYY - hXY*hXY)/det;
    
    // offset
    ofstX = -aa*gX-ab*gY-ac*gS; // in position
    ofstY = -ab*gX-bb*gY-bc*gS;
    ofstS = -ac*gX-bc*gY-cc*gS;
    /** Compute the DoG value offset */
    ofstVal = + 0.5*(gX*ofstX+gY*ofstY+gS*ofstS); //... and value

    /** output */
    *di = ofstX;
    *dj = ofstY;
    *ds = ofstS;
    *val = stack[s*w*h+i*w+j] + ofstVal;
}




/** @brief Refine the position of candidate keypoints
 *         
 *  in  @param dog : difference of Gaussian
 *  in  @param keys : list candidate keypoints (3d discrete extrema)
 * 
 *  out @param keysInterpol : list of interpolated keypoints.
 *  out @param keysReject   : list of removed keypoints.
 * 
 *  The interpolation model consists in a local 3d quadratic model of the DoG space.
 *  
 *  An interpolation is successful if the interpolated extremum lays inside the pixel area
 *  
 *  Iterative process 
 *
 */
// static void keypoints_interpolate_position(struct sift_scalespace *d, // TODO (lib_dense_anatomy)
void keypoints_interpolate_position(struct sift_scalespace *d,
                                    struct sift_keypoints *keys,
                                    struct sift_keypoints *keysInterpol,
                                    int itermax)
{

    int nMaxIntrp = itermax; /* Maximum number of consecutive unsuccessful interpolation */
    //_myfloat ofstMax = 0.6;
    _myfloat ofstMax = 0.5;
    // Ratio between two consecutive scales in the scalespace
    // assuming the ratio is constant over all scales and over all octaves
    _myfloat sigmaratio = d->octaves[0]->sigmas[1]/d->octaves[0]->sigmas[0];

    for(int k = 0; k<keys->size; k++){

        /* Loading keypoint and associated octave */
        struct keypoint* key = keys->list[k];
        int o = key->o;
        int s = key->s;
        int i = key->i;
        int j = key->j;
        struct octa* octave = d->octaves[o];
        int w = octave->w;
        int h = octave->h;
        int ns = octave->nSca;    // WARNING this includes the auxiliary scales.
        _myfloat delta  = octave->delta;
        _myfloat* imStck = octave->imStack;
        _myfloat val = key->val;

        int ic=i;   /* current value of i coordinate - at each interpolation */
        int jc=j;
        int sc=s;
        int nIntrp = 0;
        bool isConv = false;
        _myfloat ofstX = 0.;
        _myfloat ofstY = 0.;
        _myfloat ofstS = 0.;


        // ADDED 16 decembre 2014 / uy
        if (nMaxIntrp == 0){ // no interpolation
            isConv = true;
        }
        // END ADDED

        while( nIntrp < nMaxIntrp ){

            nIntrp += 1;

            /** Extrema interpolation via a quadratic function */
            /*   only if the detection is not too close to the border (so the discrete 3D Hessian is well defined) */
            if((0 < ic)&&( ic < (h-1))&&(0 < jc)&&(jc < (w-1))){
                inverse_3D_Taylor_second_order_expansion(imStck, w, h, ns, ic, jc, sc, &ofstX, &ofstY, &ofstS, &val);
            }else{
                isConv = false;
                //debug("trying interpolation on the edge of the image: iter=%i", nIntrp); // It should not happen on the first iteration
                nIntrp = nMaxIntrp;  // If the point is on the edge stop interpolation process.
                break;
            }
            /** Test if the quadratic model is consistent */
            if( (ABS(ofstX) < ofstMax) && (ABS(ofstY) < ofstMax) && (ABS(ofstS) < ofstMax) ){
                isConv = true;
                break;
            }else{ // move to another point
                // space...
                if((ofstX > +ofstMax) && ((ic+1) < (h-1))) {ic +=1;}
                if((ofstX < -ofstMax) && ((ic-1) >  0   )) {ic -=1;}
                if((ofstY > +ofstMax) && ((jc+1) < (w-1))) {jc +=1;}
                if((ofstY < -ofstMax) && ((jc-1) >  0   )) {jc -=1;}
                // ... and scale.
                if((ofstS > +ofstMax) && ((sc+1) < (ns-1))) {sc +=1;} // WARNING This must be commented for GRADUAL binary
                if((ofstS < -ofstMax) && ((sc-1) >    0  )) {sc -=1;}
            }
        }

        if(isConv == true){
            /** Create key and save in corresponding keypoint structure */
            struct keypoint* keycp = sift_malloc_keypoint_from_model_and_copy(key);
            keycp->x = (ic+ofstX)*delta;
            keycp->y = (jc+ofstY)*delta;
            keycp->i = ic;
            keycp->j = jc;
            keycp->s = sc;
            keycp->sigma = octave->sigmas[sc]*pow(sigmaratio,ofstS); /* logarithmic scale */
            keycp->val = val;
            // EXTRA - DENSE number of iteration
            keycp->iter = nIntrp;
            // EXTRA - offset for last interpolation
            keycp->ofstX = ofstX;
            keycp->ofstY = ofstY;
            keycp->ofstS = ofstS;
            // save
            sift_add_keypoint_to_list(keycp,keysInterpol);
        }
    }
}


void keypoints_interpolate_position_controlled(struct sift_scalespace *d,
                                    struct sift_keypoints *keys,
                                    struct sift_keypoints *keysInterpol,
                                    int itermax,
                                    _myfloat ofstMax_X,
                                    _myfloat ofstMax_S,
                                    int flag_jumpinscale,
                                    _myfloat fnspo,
                                    _myfloat sigma_min)
{

    int nMaxIntrp = itermax; /* Maximum number of consecutive unsuccessful interpolation */
    // Ratio between two consecutive scales in the scalespace
    // assuming the ratio is constant over all scales and over all octaves
    _myfloat sigmaratio = d->octaves[0]->sigmas[1]/d->octaves[0]->sigmas[0];

    for(int k = 0; k<keys->size; k++){

        /* Loading keypoint and associated octave */
        struct keypoint* key = keys->list[k];
        int o = key->o;
        int s = key->s;
        int i = key->i;
        int j = key->j;
        struct octa* octave = d->octaves[o];
        int w = octave->w;
        int h = octave->h;
        int ns = octave->nSca;    // WARNING this includes the auxiliary scales.
        _myfloat delta = octave->delta;
        _myfloat* imStck = octave->imStack;
        _myfloat val = key->val;

        int ic=i;   /* current value of i coordinate - at each interpolation */
        int jc=j;
        int sc=s;
        int nIntrp = 0;
        bool isConv = false;
        _myfloat ofstX = 0.;
        _myfloat ofstY = 0.;
        _myfloat ofstS = 0.;

        // Option handling to only consider the discrete extrema
        if (nMaxIntrp == 0){ // no interpolation
            isConv = true;
        }

        while( nIntrp < nMaxIntrp ){

            nIntrp += 1;

            /** Extrema interpolation via a quadratic function */
            /*   only if the detection is not too close to the border (so the discrete 3D Hessian is well defined) */
            if((0 < ic)&&( ic < (h-1))&&(0 < jc)&&(jc < (w-1))){
                inverse_3D_Taylor_second_order_expansion(imStck, w, h, ns, ic, jc, sc, &ofstX, &ofstY, &ofstS, &val);
            }else{
                isConv = false;
                nIntrp = nMaxIntrp;  // If the point is on the edge stop interpolation process.
                break;
            }
            /** Test if the quadratic model is consistent */
            if( (ABS(ofstX) < ofstMax_X) && (ABS(ofstY) < ofstMax_X) && (ABS(ofstS) < ofstMax_S) ){
                isConv = true;
                break;
            }else{ // move to another point
                // space...
                if((ofstX > +ofstMax_X) && ((ic+1) < (h-1))) {ic +=1;}
                if((ofstX < -ofstMax_X) && ((ic-1) >  0   )) {ic -=1;}
                if((ofstY > +ofstMax_X) && ((jc+1) < (w-1))) {jc +=1;}
                if((ofstY < -ofstMax_X) && ((jc-1) >  0   )) {jc -=1;}
                // ... and scale.
                if (flag_jumpinscale == 1){ 
                    if((ofstS > +ofstMax_S) && ((sc+1) < (ns-1))) {sc +=1;}
                    if((ofstS < -ofstMax_S) && ((sc-1) >    0  )) {sc -=1;}
                }
            }
        }



        // Equivalent scale of the interpolated extrema
        _myfloat scale = octave->sigmas[sc]*pow(sigmaratio,ofstS);/*logathmic scale*/

        // Check that the interpolated detection is inside the current octave 
        // (plus ofst ). This is important for non integer number of detections.
        _myfloat maxscale = sigma_min * pow(2, o + 1 + ofstMax_S/fnspo) ;
        if(( scale > maxscale )&&(nMaxIntrp > 0))
            isConv = false;
        //not the good offset for discrete. // ofstMax_S / fnspo

        if(isConv == true){
            /** Create key and save in corresponding keypoint structure */
            struct keypoint* keycp = sift_malloc_keypoint_from_model_and_copy(key);
            keycp->x = (ic+ofstX)*delta;
            keycp->y = (jc+ofstY)*delta;
            keycp->i = ic;
            keycp->j = jc;
            keycp->s = sc;
            keycp->sigma = scale; 
            keycp->val = val;
            // EXTRA - DENSE number of iteration
            keycp->iter = nIntrp;
            // EXTRA - offset for last interpolation
            keycp->ofstX = ofstX;
            keycp->ofstY = ofstY;
            keycp->ofstS = ofstS;
            sift_add_keypoint_to_list(keycp,keysInterpol);
        }
    }
}










































/** @brief  Compute Edge response 
 *    i.e.  Compute the ratio of principal curvatures
 *    i.e.  Compute the ratio (hXX + hYY)*(hXX + hYY)/(hXX*hYY - hXY*hXY);
 *          
 * 
 *     The 2D hessian of the DoG operator is computed via finite difference schemes.
 * 
 *    in  @param dog  : difference of Gaussians
 *    out @param keys : candidate keypoints
 * 
 * Note:
 *  - No keypoint is removed here
 *
 */
//STATIC
void keypoints_compute_edge_response(struct sift_scalespace *d, struct sift_keypoints *keys)
{
    for(int k=0;k<keys->size;k++){
        /* Loading keypoint and associated octave */
        struct keypoint* key = keys->list[k];
        int o = key->o;
        int s = key->s;
        int i = key->i;
        int j = key->j;
        struct octa* octave = d->octaves[o];
        int w = octave->w;
        int h = octave->h;
        _myfloat* im = &octave->imStack[s*w*h];
        /* Compute the 2d Hessian at pixel (i,j) */
        _myfloat hXX = im[(i-1)*w+j]+im[(i+1)*w+j]-2*im[i*w+j];
        _myfloat hYY = im[i*w+(j+1)]+im[i*w+(j-1)]-2*im[i*w+j] ;
        _myfloat hXY = 1./4*((im[(i+1)*w+(j+1)]-im[(i+1)*w+(j-1)])-(im[(i-1)*w+(j+1)]-im[(i-1)*w+(j-1)]));
        /* Harris and Stephen Edge response */
        _myfloat edgeResp = (hXX + hYY)*(hXX + hYY)/(hXX*hYY - hXY*hXY);
        key->edgeResp = edgeResp; /* Harris and Stephen computed on the DoG operator */
    }
}




/** @brief Dissize keys with edge response
 * 
 *  in  @param keysIn : list of keypoints, the edge response is stored
 * 
 *  out @param keysAccepted : passing
 *  out @param keysRejected : failing
 * 
 *  @param threshold on (hXX + hYY)*(hXX + hYY)/(hXX*hYY - hXY*hXY)
 *                   on the ratio of principal curvatures.
 * 
 * 
 */
//static void keypoints_discard_on_edge(struct sift_keypoints *keysIn, // TODO (lib_dense_anatomy)
void keypoints_discard_on_edge(struct sift_keypoints *keysIn,
                                      struct sift_keypoints *keysAccept,
                                      _myfloat thresh)
{
    for( int k = 0; k < keysIn->size; k++){
        struct keypoint *key = keysIn->list[k];
        bool isAccepted = ( ABS(key->edgeResp) <=  thresh);
        if (isAccepted == true){
            struct keypoint *copy = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy, keysAccept);
        }
    }
}




/** @brief Attribute a reference orientation to each keypoint of a list
 * 
 *  
 *   in  @param dx_scalespace, x component (up-bottom)
 *   in  @param dy_scalespace, y component (left-right) 
 * 
 *   in  @param keysIn     list of keypoints (pointer)
 * 
 *   out @param keysOut   list of oriented keypoints
 *
 *            size(keysIn) <= size(keysOut)          
 *
 *   @param t  : threshold over a local maxima to constitute a reference orientation 
 * 
 *   @param lambda_ori : size parameter for the Gaussian window
 * 
 * 
 * 
 */
//static void keypoints_attribute_orientations(const struct sift_scalespace *sx, // TODO (lib_dense_anatomy)
void keypoints_attribute_orientations(const struct sift_scalespace *sx,
                                             const struct sift_scalespace *sy,
                                             const struct sift_keypoints *keysIn,
                                             struct sift_keypoints *keysOut,
                                             int n_bins, _myfloat lambda_ori, _myfloat t)
{
    for(int k=0;k<keysIn->size;k++){

        // load keypoint coordinates
        struct keypoint* key = keysIn->list[k];
        _myfloat x = key->x;
        _myfloat y = key->y;
        _myfloat sigma = key->sigma;
        int o = key->o;
        int s = key->s;

        // load scalespace gradient
        int w = sx->octaves[o]->w;
        int h = sx->octaves[o]->h;
        _myfloat delta  = sx->octaves[o]->delta;
        const _myfloat* dx = &(sx->octaves[o]->imStack[s*w*h]);
        const _myfloat* dy = &(sy->octaves[o]->imStack[s*w*h]);

        //conversion to the octave's coordinates
        x /= delta;
        y /= delta;
        sigma /= delta;

        /** Accumulate gradient orientation histogram */
        sift_accumulate_orientation_histogram(x, y, sigma, dx, dy, w, h, n_bins, lambda_ori, key->orihist);

        /** Extract principal orientation */
        _myfloat* principal_orientations = xmalloc(n_bins*sizeof(_myfloat));
        int n_prOri;
        n_prOri = sift_extract_principal_orientations(key->orihist, n_bins, t, principal_orientations); /*t = 0.8 threhsold for secondary orientation */

        /** Updating keypoints and save in new list */
        for(int n = 0; n < n_prOri; n++){
            struct keypoint* copy = sift_malloc_keypoint_from_model_and_copy(key);
            copy->theta = principal_orientations[n];
            sift_add_keypoint_to_list(copy, keysOut);
        }
        free(principal_orientations);
    }
}




static void keypoints_attribute_one_orientation(const struct sift_scalespace *sx,
                                                const struct sift_scalespace *sy,
                                                struct sift_keypoints *keys,
                                                int n_bins, _myfloat lambda_ori)
{
    for(int k = 0; k < keys->size; k++){

        // load keypoint coordinates
        struct keypoint* key = keys->list[k];
        _myfloat x = key->x;
        _myfloat y = key->y;
        _myfloat sigma = key->sigma;
        int o = key->o;
        int s = key->s;

        // load scalespace gradient
        int w = sx->octaves[o]->w;
        int h = sx->octaves[o]->h;
        _myfloat delta  = sx->octaves[o]->delta;
        const _myfloat* dx = &(sx->octaves[o]->imStack[s*w*h]);
        const _myfloat* dy = &(sy->octaves[o]->imStack[s*w*h]);

        //conversion to the octave's coordinates
        x /= delta;
        y /= delta;
        sigma /= delta;

        /** Accumulate gradient orientation histogram */
        sift_accumulate_orientation_histogram(x, y, sigma, dx, dy, w, h, n_bins,lambda_ori,key->orihist);

        /** Extract principal orientation (includes histogram smoothing)*/
        key->theta = sift_extract_one_orientation(key->orihist, key->n_bins);
    }
}



//static void keypoints_discard_near_the_border(struct sift_keypoints *keysIn,
void keypoints_discard_near_the_border(struct sift_keypoints *keysIn,
                                            struct sift_keypoints *keysAccept,
                                            int w,
                                            int h,
                                            _myfloat lambda)
{
    for( int k = 0; k < keysIn->size; k++){
        struct keypoint *key = keysIn->list[k];
        _myfloat x = key->x;
        _myfloat y = key->y;
        _myfloat sigma = key->sigma;
        bool isAccepted = (x-lambda*sigma > 0.0 )&&( x + lambda*sigma < (_myfloat)h)
                       && (y-lambda*sigma > 0.0 )&&( y + lambda*sigma < (_myfloat)w);
        if (isAccepted == true){
            struct keypoint *copy = sift_malloc_keypoint_from_model_and_copy(key);
            sift_add_keypoint_to_list(copy, keysAccept);
        }
    }
}




/** @brief Attribute a feature vector to each keypoint of a list
 *
 *
 *   in  @param dx_scalespace, x component (up to bottom)
 *   in  @param dy_scalespace, y component (left to right)
 *
 *   in  @param keys     list of keypoints to be described
 *
 *   @param n_hist   , the descriptor is an array of nhist^2 weighted histograms.
 *   @param n_ori    , each weighted histogram has n_ori bins.
 *
 *   @param lambda_descr : size parameter for the Gaussian window
 *                         - the Gaussian window has a std dev of lambda_descr X sigma
 *                         - the patch considered has a w of
 *                         2*(1+1/n_hist)*lambda_descr X sigma
 */
//static void keypoints_attribute_descriptors(struct sift_scalespace *sx,
void keypoints_attribute_descriptors(struct sift_scalespace *sx,
                                            struct sift_scalespace *sy,
                                            struct sift_keypoints *keys,
                                            int n_hist,
                                            int n_ori,
                                            _myfloat lambda_descr)
{
    int n_descr = n_hist*n_hist*n_ori;

    for(int k = 0; k < keys->size; k++){

        /** Loading keypoint gradient scalespaces */
        struct keypoint* key = keys->list[k];
        _myfloat x = key->x;
        _myfloat y = key->y;
        int o = key->o;
        int s = key->s;
        _myfloat sigma = key->sigma;
        _myfloat theta = key->theta;

        // load scalespace gradient
        int w = sx->octaves[o]->w;
        int h = sx->octaves[o]->h;
        _myfloat delta  = sx->octaves[o]->delta;
        const _myfloat* dx = &(sx->octaves[o]->imStack[s*w*h]);
        const _myfloat* dy = &(sy->octaves[o]->imStack[s*w*h]);

        // conversion to the octave's coordinates
        x /= delta;
        y /= delta;
        sigma /= delta;

        /** Compute descriptor representation */
        sift_extract_feature_vector(x, y, sigma, theta,
                                    dx, dy, w, h,
                                    n_hist, n_ori, lambda_descr,
                                    key->descr);

        /** Threshold and quantization of the descriptor */
        sift_threshold_and_quantize_feature_vector(key->descr, n_descr, 0.2);
    }
}





struct sift_parameters* sift_assign_default_parameters()
{
    struct sift_parameters* p = xmalloc(sizeof(*p));
    p->n_oct = 8;
    p->n_spo = 3;
    p->sigma_min = 0.8;
    p->delta_min = 0.5;
    p->sigma_in = 0.5;
    p->C_DoG = 0.013333333;  // = 0.04/3
    p->C_edge = 10;
    p->n_bins = 36;
    p->lambda_ori = 1.5;
    p->t = 0.80;
    p->n_hist = 4;
    p->n_ori = 8;
    p->lambda_descr = 6;
    p->itermax = 5;

    // EXTRA
    p->dog_nspo = 3; // definition of the DoG operator.
    p->epsilon = FLT_EPSILON; // definition of a discrete extrema
    p->fnspo = 3.0;  // density of scales per octave (non integer !)
    // interpolation process.
 //   p->ofstMax_X = 0.5;
 //   p->ofstMax_S = 0.5;
    p->ofstMax_X = 0.6;
    p->ofstMax_S = 0.6;

    p->flag_jumpinscale = 1;
    // controlling extraction of 3d extrema
    p->discrete_extrema_dR = 1;
    p->discrete_sphere_r = 1;
    p->discrete_sphere_dr = 1;
    // to confirm that there is an extremum
    p->ball_r1 = 1;
    p->ball_r2 = 2;
    p->ball_r3 = 3;
    p->ball_alpha = 0.24;
    p->ball_alpha = 0.8;

    return p;
}


// The number of octaves is limited by the size of the input image
// STATIC
int number_of_octaves(int w, int h, const struct sift_parameters* p)
{
    // The minimal size (width or height) of images in the last octave.
    int hmin = 12;
    // The size (min of width and height) of images in the first octave.
    int h0 = MIN(w,h)/p->delta_min;
    // The number of octaves.
    int n_oct = MIN(p->n_oct, (int)(log(h0/hmin)/M_LN2) + 1);
    return n_oct;
}


// To make the threshold independent of the scalespace discretization along
// scale (number of scale per octave).
//STATIC
_myfloat convert_threshold(const struct sift_parameters* p)
{
    // converting the threshold to make it consistent
    //_myfloat k_nspo =  exp( M_LN2/( (_myfloat)(p->n_spo)));
    _myfloat k_nspo =  exp( M_LN2/( (_myfloat)(p->n_spo)));
    _myfloat k_3 =  exp( M_LN2/( (_myfloat)3));
    _myfloat thresh = (k_nspo - 1) / (k_3 - 1) * p->C_DoG ;
    return thresh;
}



struct sift_keypoints* sift_anatomy(const _myfloat* x, int w, int h, const struct sift_parameters* p,
                                    struct sift_scalespace* ss[4],
                                    struct sift_keypoints* kk[6])

{
    // WARNING
    //  ss[2]: pointers to two scalespace structures for lates use (The Gaussian
    //         scale-space and the DoG scale-space).
    //  kk[6]: pointers to six keypoint lists structure for later use.
    
    struct sift_keypoints* k = sift_malloc_keypoints();
   
    // The number of octaves.
    int n_oct = number_of_octaves(w, h, p);
    
    // adapt the threshold to the scalespace discretization
    _myfloat thresh = convert_threshold(p);
    
    /** MEMORY ALLOCATION **/
    /** scale-space structure */
    struct sift_scalespace* s   = sift_malloc_scalespace_lowe(n_oct,p->n_spo,w,h,p->delta_min,p->sigma_min);
    struct sift_scalespace* d   = sift_malloc_scalespace_dog_from_scalespace(s);
    struct sift_scalespace* sx  = sift_malloc_scalespace_from_model(s);
    struct sift_scalespace* sy  = sift_malloc_scalespace_from_model(s);
    /** list-of-keypoints (Already allocated) */
    struct sift_keypoints* kA   = kk[0];  /* 3D (discrete) extrema   */
    struct sift_keypoints* kB   = kk[1];  /* passing the threshold on DoG  */
    struct sift_keypoints* kC   = kk[2];  /* interpolated 3D extrema (continuous) */
    struct sift_keypoints* kD   = kk[3];  /* passing the threshold on DoG  */
    struct sift_keypoints* kE   = kk[4];  /* passing OnEdge filter */
    struct sift_keypoints* kF   = kk[5];  /* keypoints whose distance to image borders are sufficient */

    /** KEYPOINT DETECTION ***************************************************/
    scalespace_compute(s, x, w, h, p->sigma_in);  /* Builds Lowe's scale-space */
    scalespace_compute_dog(s,d);

    keypoints_find_3d_discrete_extrema(d, kA, p->n_ori, p->n_hist, p->n_bins);
    keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
    keypoints_interpolate_position(d, kB, kC, p->itermax);
    //keypoints_interpolate_positionTEMP(d, kB, kC, p->itermax, p->ofstMax_X, p->ofstMax_S, p->flag_jumpinscale);
    keypoints_discard_with_low_response(kC, kD, thresh); 
    keypoints_compute_edge_response(d,kD); 
    keypoints_discard_on_edge(kD, kE, (p->C_edge+1)*(p->C_edge+1)/p->C_edge);
    keypoints_discard_near_the_border(kE, kF, w, h, 1.0);


    /** KEYPOINT DESCRIPTION *************************************************/
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



struct sift_keypoints* sift_anatomy_without_description(const _myfloat* x, int w, int h,
                                                        const struct sift_parameters* p,
                                                        struct sift_scalespace* ss[2], 
                                                        struct sift_keypoints* kk[5])


{
    // WARNING
    //  ss[2]: pointers to two scalespace structures for lates use (The Gaussian
    //         scale-space and the DoG scale-space).
    //  kk[5]: pointers to five keypoint lists structure for later use.

    struct sift_keypoints* k = sift_malloc_keypoints();

    // The number of octaves.
    int n_oct = number_of_octaves(w, h, p);
    
    // Adapt the threshold to the scalespace discretization
    _myfloat thresh = convert_threshold(p);

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

    /** KEYPOINT DETECTION ***************************************************/
    scalespace_compute(s, x, w, h, p->sigma_in);  /* Builds Lowe's scale-space */
    scalespace_compute_dog(s,d);
    keypoints_find_3d_discrete_extrema(d, kA, p->n_ori, p->n_hist, p->n_bins);
    keypoints_discard_with_low_response(kA, kB, 0.8*thresh);
    keypoints_interpolate_position(d, kB, kC, p->itermax);
    keypoints_discard_with_low_response(kC, kD, thresh);
    keypoints_compute_edge_response(d,kD);
    keypoints_discard_on_edge(kD, kE, (p->C_edge+1)*(p->C_edge+1)/p->C_edge);
    keypoints_discard_near_the_border(kE, k, w, h, 1.0);

    /** scalespace structures*/
    ss[0] = s;
    ss[1] = d;

    return k;
}




void sift_anatomy_only_description(const _myfloat* x, int w, int h, const struct sift_parameters* p,
                                    struct sift_keypoints* k)
{
    // The number of octaves.
    int n_oct = number_of_octaves(w, h, p);

    // MEMORY ALLOCATION
    struct sift_scalespace* s   = sift_malloc_scalespace_lowe(n_oct,p->n_spo,w,h,p->delta_min,p->sigma_min);
    struct sift_scalespace* sx  = sift_malloc_scalespace_from_model(s);
    struct sift_scalespace* sy  = sift_malloc_scalespace_from_model(s);

    scalespace_compute(s, x, w, h, p->sigma_in);  /* Builds Lowe's scale-space */
    scalespace_compute_gradient(s,sx,sy); /* Pre-computes gradient scale-space */
    keypoints_attribute_descriptors(sx, sy, k, p->n_hist,p->n_ori,p->lambda_descr);

    sift_free_scalespace(s);
    sift_free_scalespace(sx);
    sift_free_scalespace(sy);
}







void sift_anatomy_orientation_and_description(const _myfloat* x,
                                              int w,
                                              int h,
                                              const struct sift_parameters* p,
                                              struct sift_keypoints* k)
{
    // The number of octaves.
    int n_oct = number_of_octaves(w, h, p);

    // MEMORY ALLOCATION
    struct sift_scalespace* s   = sift_malloc_scalespace_lowe(n_oct,p->n_spo,w,h,p->delta_min,p->sigma_min);
    struct sift_scalespace* sx  = sift_malloc_scalespace_from_model(s);
    struct sift_scalespace* sy  = sift_malloc_scalespace_from_model(s);
    scalespace_compute(s, x, w, h, p->sigma_in);  // Builds Lowe's scale-space
    scalespace_compute_gradient(s,sx,sy); // Pre-computes gradient scale-space

    keypoints_attribute_one_orientation(sx, sy, k, p->n_bins, p->lambda_ori);
    keypoints_attribute_descriptors(sx, sy, k, p->n_hist,p->n_ori,p->lambda_descr);

    sift_free_scalespace(s);
    sift_free_scalespace(sx);
    sift_free_scalespace(sy);
}
