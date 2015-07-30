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

    [3] "Method and apparatus for identifying scale invariant features
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


This program is free software: you can use, modify and/or
redistribute it under the terms of the simplified BSD
License. You should have received a copy of this license along
this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.

*/
/**
 * @file sift_description.c
 * @brief Computation the SIFT feature vector
 *
 * @li Attribution of a principal orientation
 * @li Computation of the SIFT feature vector
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */

#ifndef _LIB_DESCRIPTION_H_
#define _LIB_DESCRIPTION_H_

#include "lib_util.h"  // for the quadruple precision float

void sift_accumulate_orientation_histogram(_myfloat x_key,
                                           _myfloat y_key,
                                           _myfloat sigma_key,
                                           const _myfloat* gradX,
                                           const _myfloat* gradY,
                                           int w,
                                           int h,
                                           int nbins,
                                           _myfloat lambda_ori,
                                           _myfloat* hist);

// Apply smoothing and extract local extrema
int sift_extract_principal_orientations(_myfloat* hist,
                                        int nbins,
                                        _myfloat threshold,
                                        _myfloat* principal_orientations);

// Apply smoothing and extract the global extremum
_myfloat sift_extract_one_orientation(_myfloat* hist, int nbins);

// Accumulate the descriptor's histograms
void sift_extract_feature_vector(_myfloat x_key,
                                _myfloat y_key,
                                _myfloat sigma_key,
                                _myfloat theta_key,
                                const _myfloat* gradX,
                                const _myfloat* gradY,
                                int w,
                                int h,
                                int Nhist,
                                int Nbins,
                                _myfloat lambda_descr,
                                _myfloat* descr);

void sift_threshold_and_quantize_feature_vector(_myfloat* descr, int n, _myfloat threshold);

#endif // _LIB_DESCRIPTION_H_
