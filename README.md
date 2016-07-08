# sift_anatomy_extra

## Overview

Implementation of the SIFT method.
This C/C++ ANSI source code is related to the article

    1. "An Analysis of the scale-space sampling in SIFT."
        I. Rey Otero, J.M. Morel,  M. Delbracio
        ICIP 2014.
        http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=7025982&tag=1

Part of the code is also related to the article

    2. "Anatomy of the SIFT method."
        I. Rey Otero, M. Delbracio
        IPOL (Image Processing On Line - 2014)
        http://www.ipol.im/pub/art/2014/82/
        http://dx.doi.org/10.5201/ipol.2014.82


This adds extra functionalities for the numerical analysis of the SIFT method.
 - Exact scale-space computation using DFT interpolation
 - spline interpolation of the input image
 - extrema interpolation satisfying the least square error
 - EXR input/output
 - define the DoG differential operator independently of the scale-space
  sampling

## Patent Warning and License

The methods implemented here are extensions of the SIFT method.
This method is patented

    3. "Method and apparatus for identifying scale invariant features in an image."
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

## Compiling (Linux)

To compile, type
```
$ make
```
in the directory where the Makefile is located.
The code relies of the following common libraries:
  - libpng12-dev
  - libjpeg-dev
  - libtiff5-dev
  - libfftw3-dev
  - libopenexr-dev

The compilation of the source code provides three executables:

1. sift_cli:  (original algorithm) applies the SIFT method to a PNG image.
              Parameters are documented in [2]. 

2. match_cli:   matches the SIFT keypoints extracted from two images.

3. extra_sift_cli:  Compute the Gaussian scale-space using an exact
                    implementation of the Gaussian convolution based on the
                    Fourier interpolation.
                    The DoG operator is defined independently of the 
                    scale-space sampling.

4. gradual: same as extra_sift_cli, computes the scale-space gradually to
            reduce the memory usage.


## SIFT executable - Usage 

```
$ ./sift_cli image [options...]  [> keys]
```

options :

    -ss_noct        (8)  number of octaves
    -ss_nspo        (3)  number of scales per octaves
    -ss_dmin      (0.5)  the sampling distance in the first octave
    -ss_smin      (0.8)  blur level on the seed image
    -ss_sin       (0.5)  assumed level of blur in the input image

    -thresh_dog (0.0133) threshold over the DoG response
    -thresh_edge   (10)  threshold over the ratio of principal curvature

    -ori_nbins    (36)   number of bins in the orientation histogram
    -ori_thresh  (0.8)   threhsold for considering local maxima in
                         the orientation histogram
    -ori_lambda  (1.5)   sets how local is the analysis of the gradient
                         distribution

    -descr_nhist   (4)   number of histograms per dimension
    -descr_nori    (8)   number of bins in each histogram
    -descr_lambda  (6)   sets how local the descriptor is

    -verb_keys   label   flag to output the intermediary sets of keypoints
    -verb_ss     label   flag to output the scalespaces (Gaussian and DoG)

EXTRA parameters

    -ss_fnspo (3.00) float number of scales per octaves (-ss_nspo not used) 
    -flag_semigroup   BOOL (1)    semigroup (1) or direct (0)               
    -flag_dct         BOOL (0)    dct (1) or discrete (0)                   
    -flag_log         BOOL (0)    normalized Laplacian (1) or DoG (0)       
    -flag_interp     (0)  bilin (0) / DCT (1)/ bsplines (3,5,..,11)             
    -itermax             5        max number of iterations                  
    -epsilon        FLT_EPSILON    (for _myfloat comparison)                
    -dog_nspo    3  number of scales per octaves for DoG operator definition
    -ofstMax_X   (0.5)   interpolation validity domain definition in space  
    -ofstMax_S   (0.5)                 ... in scale                         
    -flag_jumpinscale   (0 in gradual / not an option in gradual)  


