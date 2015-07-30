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
 * @file lib_sift_anatomy.h
 * @brief SIFT anatomy interface
 *
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */


#ifndef _LIB_SIFT_ANATOMY_H_
#define _LIB_SIFT_ANATOMY_H_


#include "lib_scalespace.h"
#include "lib_keypoint.h"
#include "lib_util.h"


struct sift_parameters
{
    int n_oct, n_spo, n_hist, n_bins, n_ori, itermax;
    _myfloat sigma_min, delta_min, sigma_in, C_DoG, C_edge, lambda_ori, t, lambda_descr;
    // EXTRA
    int dog_nspo;
    _myfloat epsilon;
    _myfloat fnspo;  // float value  for the number of scales per octave (only gradual).
    _myfloat ofstMax_X; // for extrema interpolation (in space).
    _myfloat ofstMax_S; // for extrema interpolation (in the scale direction).
    int flag_jumpinscale;
};






struct sift_parameters* sift_assign_default_parameters();

struct sift_keypoints* sift_anatomy(const _myfloat* x, int w, int h, const struct sift_parameters* p,
                                    struct sift_scalespace* ss[4],
                                    struct sift_keypoints* kk[6]);

struct sift_keypoints* sift_anatomy_without_description(const _myfloat* x, int w, int h, const struct sift_parameters* p,
                                    struct sift_scalespace* ss[2],
                                    struct sift_keypoints* kk[5]);

void sift_anatomy_only_description(const _myfloat* x, int w, int h, const struct sift_parameters* p, struct sift_keypoints* k);

void sift_anatomy_orientation_and_description(const _myfloat* x, int w, int h, const struct sift_parameters* p, struct sift_keypoints* k);


void scalespace_compute(struct sift_scalespace* ss,
                               const _myfloat* image,
                               int im_w,
                               int im_h,
                               _myfloat sigma_in);

void scalespace_compute_gradient(const struct sift_scalespace* scalespace,
                                 struct sift_scalespace* sx,
                                 struct sift_scalespace* sy);






// TODO TEMP (These are added for lib_dense_anatomy.c, they are statically declared otherwise)

int number_of_octaves(int w, int h, const struct sift_parameters* p);

_myfloat convert_threshold(const struct sift_parameters* p);


void scalespace_compute_dog(const struct sift_scalespace *s,
                                   struct sift_scalespace *d);

void keypoints_find_3d_discrete_extrema(struct sift_scalespace* d,
                                               struct sift_keypoints* keys,
                                               int n_ori,
                                               int n_hist,
                                               int n_bins);

void keypoints_find_3d_discrete_extrema_epsilon(struct sift_scalespace* d,
                                               struct sift_keypoints* keys,
                                               int n_ori,
                                               int n_hist,
                                               int n_bins,
                                               _myfloat epsilon);


void keypoints_interpolate_position(struct sift_scalespace *d,
                                    struct sift_keypoints *keys,
                                    struct sift_keypoints *keysInterpol,
                                    int itermax);




void keypoints_discard_with_low_response(struct sift_keypoints *keysIn,
                                                struct sift_keypoints *keysAccept,
                                                _myfloat thresh);

void keypoints_attribute_orientations(const struct sift_scalespace *sx,
                                             const struct sift_scalespace *sy,
                                             const struct sift_keypoints *keysIn,
                                             struct sift_keypoints *keysOut,
                                             int n_bins, _myfloat lambda_ori, _myfloat t);

void keypoints_discard_near_the_border(struct sift_keypoints *keysIn,
                                            struct sift_keypoints *keysAccept,
                                            int w,
                                            int h,
                                            _myfloat lambda);

void keypoints_compute_edge_response(struct sift_scalespace *d,
        struct sift_keypoints *keys);


void keypoints_discard_on_edge(struct sift_keypoints *keysIn,
                                      struct sift_keypoints *keysAccept,
                                      _myfloat thresh);

void keypoints_attribute_descriptors(struct sift_scalespace *sx,
                                            struct sift_scalespace *sy,
                                            struct sift_keypoints *keys,
                                            int n_hist,
                                            int n_ori,
                                            _myfloat lambda_descr);


// TO control the interpolation iterative process and avoid border effects in the case of non integer nspo
void keypoints_interpolate_position_controlled(struct sift_scalespace *d,
                                        struct sift_keypoints *keys,
                                        struct sift_keypoints *keysInterpol,
                                        int itermax,
                                        _myfloat ofstMax_X,
                                        _myfloat ofstMax_S,
                                        int flag_jumpinscale,
                                        _myfloat fnspo,
                                        _myfloat sigma_min);

#endif // _LIB_SIFT_ANATOMY_H_
