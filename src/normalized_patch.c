#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <float.h>
#include <stdarg.h> // for debug message

#include "lib_sift_anatomy.h"
#include "lib_ellipse_anatomy.h"
#include "lib_util.h"
#include "io_png.h"

#define MAX(i,j) ( ((i)<(j)) ? (j):(i) )
#define MIN(i,j) ( ((i)<(j)) ? (i):(j) )
#define ABS(i)   ( ((i)>0 ) ? (i):(-(i)) )


_myfloat* read_image_to_gray(char* fname, int *pw, int *ph)
{
    size_t w, h;
    _myfloat* imf32 = io_png_read_f32_rgb(fname, &w, &h);
    _myfloat* x = xmalloc(w*h*sizeof(_myfloat));
    for(unsigned int i = 0; i < w*h; i++)
        x[i] = (_myfloat)((0.299*imf32[i] + 0.587*imf32[w*h+i] + 0.1140*imf32[2*w*h+i])/256.); //RGB2GRAY
    xfree(imf32);
    *pw = w;
    *ph = h;
    return x;
}


void print_usage()
{
    fprintf(stderr, "SIFT with affine normalized keypoints ver 20140823         \n");
    fprintf(stderr, "Usage:  normalized_patch ellipsesfile image [lambda_descr gaus_flag] \n");
}

int main(int argc, char **argv)
{

    if((argc != 3)&&(argc != 5)){
        print_usage();
        return EXIT_FAILURE;
    }

    // Loading image
    int w, h;
    _myfloat* im = read_image_to_gray(argv[2], &w, &h);

    // load standard method parameter
    struct sift_parameters* p = sift_assign_default_parameters();
    bool gauss_flag = true;
    if (argc  == 5){
        p->lambda_descr = atof(argv[argc-2]);
        gauss_flag = atoi(argv[argc-1]);
    }

    // load keypoints (ellipse)
    struct sift_keypoints *l = sift_malloc_keypoints();
    sift_read_ellipses(l, argv[1], p->n_hist, p->n_ori, p->n_bins, 0);
    fprintf(stderr," after read");
    estimate_scalespace_index(l, p);
    fprintf(stderr," after estimate scalespace index");

    // scalespace, gradient, orientation and description
    ellipse_sift_anatomy_orientation_and_description(im, w, h, p, l, gauss_flag);
    fprintf(stderr," after ellipse index");

    // printing ellipses
    fprintf_ellipses(stdout, l, 1);

    // memory deallocation
    sift_free_keypoints(l);
    free(im);

    return EXIT_SUCCESS;
}
