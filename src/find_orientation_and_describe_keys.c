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





== Patent Warning and Licence =================================================

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
 * @file sift.c
 * @brief [[MAIN]] The SIFT method 
 *
 * @li basic SIFT transform applied to one image
 * @li verbose SIFT transform 
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "lib_sift.h"
#include "lib_sift_anatomy.h"
#include "lib_io_scalespace.h"
#include "io_png.h"




static void *xmalloc(size_t size)
{
    if (size == 0)
        fprintf(stderr,"xmalloc: zero size");
        void *p = malloc(size);
        if (!p)
        {
            double sm = size / (0x100000 * 1.0);
            fprintf(stderr,"xmalloc: out of memory when requesting "
            "%zu bytes (%gMB)",//:\"%s\"",
            size, sm);//, strerror(errno));
        }
        return p;
}


static void xfree(void *p)
{
    if (!p)
        fprintf(stderr,"trying to free a null pointer");
        free(p);
}



void print_usage()
{
    fprintf(stderr, "HELP:  compute the SIFT descriptors on keypoints provided by the user\n\n");
    fprintf(stderr, "     describe_keys keys image \n");
    fprintf(stderr, "               output: x y scale orientation descriptor[128]   in standard output\n\n");

}




float* read_image_to_gray(char* fname, int *pw, int *ph)
{
    size_t w, h;
    float* imf32 = io_png_read_f32_rgb(fname, &w, &h);

    float* x = xmalloc(w*h*sizeof(float));
    for(unsigned int i = 0; i < w*h; i++)
        x[i] = (float)((0.299*imf32[i] + 0.587*imf32[w*h+i] + 0.1140*imf32[2*w*h+i])/256.); //RGB2GRAY
    xfree(imf32);
    *pw = w;
    *ph = h;
    return x;
}








/** @brief Main SIFT routine
 *
 * takes as input 1) a list of keypoints coordinates and 2) an image. 
 *
 */
int main(int argc, char **argv)
{

    if (argc !=3){
        print_usage();
        return -1;
    }

    /** Load image */
    int w, h;
    float* x = read_image_to_gray(argv[2], &w, &h);

    /** Read keypoints locations from file */ 
    int n;
    struct sift_keypoint_std * k = sift_read_keyslocation_from_file(argv[1], &n);

    /** Compute the SIFT descriptors */
    sift_find_ori_and_fill_descriptors(x, w, h, k, n);

    /** Save into a file the keypoints with their orientations and their feature vectors */
    fprintf_keypoint_std(stdout, k, n); 
    sift_write_to_file("out.tt", k, n); // TEMP

    /* memory deallocation */
    xfree(x);
    xfree(k);
    return 0;
}

