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
#include <string.h>


extern "C" {

#include "lib_sift_anatomy.h"
#include "lib_io_scalespace.h"
#include "io_png.h"
#include "lib_util.h"
#include "iio.h"
#include "lib_dense_anatomy.h"

// To read keypoints
#include "lib_sift.h"

}

#include "io_exr.h"





void print_usage()
{
    fprintf(stderr, "Anatomy EXTRACTS 33X33 PATCHES SURROUNDING EXTREMA IN A LIST    \n");
    fprintf(stderr, "Usage:  extra_sift image KEYPOINTS [options...]                          \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -ss_noct        (8)  number of octaves                                  \n");
    fprintf(stderr, "   -ss_nspo        (3)  number of scales per octaves                       \n");
    fprintf(stderr, "   -ss_dmin      (0.5)  the sampling distance in the first octave          \n");
    fprintf(stderr, "   -ss_smin      (0.8)  blur level on the seed image                       \n");
    fprintf(stderr, "   -ss_sin       (0.5)  assumed level of blur in the input image           \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -thresh_dog (0.0133) threshold over the DoG response                    \n");
    fprintf(stderr, "   -thresh_edge   (10)  threshold over the ratio of principal curvature    \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -ori_nbins    (36)   number of bins in the orientation histogram        \n");
    fprintf(stderr, "   -ori_thresh  (0.8)   threhsold for considering local maxima in          \n");
    fprintf(stderr, "                        the orientation histogram                          \n");
    fprintf(stderr, "   -ori_lambda  (1.5)   sets how local is the analysis of the gradient     \n");
    fprintf(stderr, "                        distribution                                       \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -descr_nhist   (4)   number of histograms per dimension                 \n");
    fprintf(stderr, "   -descr_nori    (8)   number of bins in each histogram                   \n");
    fprintf(stderr, "   -descr_lambda  (6)   sets how local the descriptor is                   \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -verb_keys   label   flag to output the intermediary sets of keypoints  \n");
    fprintf(stderr, "   -verb_ss     label   flag to output the scalespaces (Gaussian and DoG)  \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "          -----  EXTRA   dense_anatomy -----                               \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -flag_semigroup   BOOL (1)    semigroup (1) or direct (0)               \n");
    fprintf(stderr, "   -flag_dct         BOOL (0)    dct (1) or discrete (0)                   \n");
    fprintf(stderr, "   -flag_log         BOOL (0)    normalized Laplacian (1) or DoG (0)       \n");
    fprintf(stderr, "   -itermax             5        max number of iterations                  \n");
    fprintf(stderr, "   -flag_interp     (0)  bilin (0) / DCT (1)/ bsplines (2..11)             \n");
}


/**
 *
 * Output 
 *   -1 : malformed argument
 *    0 : option not found  
 *    1 : option found
 */
static int pick_option(int* c, char*** v, char* opt, char* val)
{
    int output = 0;
    int argc = *c;
    char **argv = *v;
    // scan the command line for '-opt'
    for(int i = 0; i < argc; i++){
        if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, opt))
        {
            // check for a corresponding value
            if (i == argc-1){
                output = -1;
            }
            else{
                if (argv[i+1][0] == '-'){
                    output  = -1;
                }
                // if the option call is well formed
                else{
                    // copy the option value ...
                    strcpy(val, argv[i+1]);
                    // ... and remove from the command line
                    for (int j = i; j < argc - 2; j++){
                        (*v)[j] = (*v)[j+2];
                    }
                    *c -= 2;
                    output = 1;
                }
            }
            // print an error if not succes
            if(output == -1){
                fprintf(stderr, "Fatal error: option %s requires an argument.\n", opt);
            }
        }
    }
    return output;
}



static int parse_options(int argc, char** argv,
                         struct sift_parameters* p,
                         int *flag_keys,
                         int *flag_ss,
                         char* label_keys,
                         char* label_ss,
                         int *flag_semigroup,
                         int *flag_dct,
                         int *flag_log,
                         int *flag_interp)
{
    int isfound;
    char val[128];

    isfound = pick_option(&argc, &argv, "ss_noct", val);
    if (isfound ==  1)    p->n_oct = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ss_nspo", val);
    if (isfound ==  1)    p->n_spo = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ss_dmin", val);
    if (isfound ==  1)    p->delta_min = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ss_smin", val);
    if (isfound ==  1)    p->sigma_min = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ss_sin", val);
    if (isfound ==  1)    p->sigma_in = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "thresh_dog", val);
    if (isfound ==  1)    p->C_DoG = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "thresh_edge", val);
    if (isfound ==  1)    p->C_edge = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ori_nbins", val);
    if (isfound ==  1)    p->n_bins = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ori_thresh", val);
    if (isfound ==  1)    p->t = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ori_lambda", val);
    if (isfound ==  1)    p->lambda_ori = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "descr_nhist", val);
    if (isfound ==  1)    p->n_hist = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "descr_nori", val);
    if (isfound ==  1)    p->n_ori = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "descr_lambda", val);
    if (isfound ==  1)    p->lambda_descr = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "verb_keys", val);
    if (isfound ==  1){
        *flag_keys = 1;
        strcpy(label_keys, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "verb_ss", val);
    if (isfound ==  1){
        *flag_ss = 1;
        strcpy(label_ss, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;

    // EXTRA DENSE ANATOMY
    isfound = pick_option(&argc, &argv, "flag_semigroup", val);
    if (isfound ==  1)    *flag_semigroup = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "flag_dct", val);
    if (isfound ==  1)    *flag_dct = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "flag_log", val);
    if (isfound ==  1)    *flag_log = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "itermax", val);
    if (isfound ==  1)    p->itermax = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "flag_interp", val);
    if (isfound ==  1)    *flag_interp = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    // check for unknown option call
    for(int i = 0; i < argc; i++){
        if (argv[i][0] == '-'){
            fprintf(stderr, "Fatal error: option \"-%s\" is unknown.\n", argv[i]+1);
            print_usage();
            return EXIT_FAILURE;
        }
    }
    // check for input image
    //if (argc != 2){
    if (argc != 3){
        print_usage();
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}




// NOTE: copied from iio.c 
static bool string_suffix(const char *s, const char *suf)
{
	int len_s = strlen(s);
	int len_suf = strlen(suf);
	if (len_s < len_suf)
		return false;
	return 0 == strcmp(suf, s + (len_s - len_suf));
}

// reads the filename suffix and uses the appropriate routine in function.
static _myfloat* read_image(char* n, int* w, int* h)
{
    float* x;

    if (string_suffix(n, ".PNG") || string_suffix(n, ".png")){
        // IO_PNG
        size_t wt, ht;
        x = io_png_read_f32_gray(n, &wt, &ht);
        for(unsigned int i=0; i < wt*ht; i++)
            x[i] /= 256.0;
        *w = wt; *h = ht;
    }
    else if ( string_suffix(n, ".EXR") || string_suffix(n, ".exr")){
        // IO_EXR
        x = readEXR_float(n, w, h);
    }
    else{
        // IIO
        x = iio_read_image_float(n, w, h);
        for(int i=0; i < (*w)*(*h); i++)
            x[i] /= 256.0;
    }

    //return x;
    
    _myfloat* y = (_myfloat*)malloc((*w)*(*h)*sizeof(_myfloat));
    for (int i = 0; i < (*w)*(*h); i++){
        y[i] = (_myfloat)x[i];
    }
    free(x);
    return(y);
}



void print_keypoints_and_vals(const struct sift_keypoints* keys, int dog_nspo)
{
    _myfloat k_nspo =  exp( M_LN2/( (_myfloat)(dog_nspo)));
    _myfloat k_3 =  exp( M_LN2/( (_myfloat)3));
    _myfloat factor = (k_3 - 1) / (k_nspo - 1);
    //
    for (int i = 0; i < keys->size; i++){
        struct keypoint* k = keys->list[i];
        _myfloat val = k->val * factor; // normalization (all relative to nspo = 3)
        fprintf(stdout, "%f %f %f %f %i \n", k->x, k->y, k->sigma, val, k->iter);
        //fprintf(stdout, "%f %f %f %f %i %i \n", k->x, k->y, k->sigma, val, k->o, k->s);
    }
}


void print_patches(_myfloat* patches, int n)
{
    for(int l = 0; l < n ; l++){
        for(int c = 0; c < 33*33; c++){
            fprintf(stdout, "%f ", patches[l*(33*33)+c]);
        }
        fprintf(stdout, "\n");
    }
}






static struct sift_keypoints* sift_translate_standard_into_anatomy_with_PARAM(const struct sift_keypoint_std* k,
                                                                       int n,
                                                                       const struct sift_parameters* p)
{
    // load the default parameters are required
    float sigma_min = p->sigma_min;  // 0.8
    float delta_min = p->delta_min;  // 0.5
    int n_spo = p->n_spo;   // 3
    int n_ori = p->n_ori;   // 8
    int n_hist = p->n_hist; // 4
    int n_bins = p->n_bins; // 36

    struct sift_keypoints* keys = sift_malloc_keypoints();
    for(int i = 0; i < n; i++){
        struct keypoint* key = sift_malloc_keypoint(n_ori, n_hist, n_bins);
        /* reading the extremum continuous coordinates */
        key->x = k[i].x;
        key->y = k[i].y;
        key->sigma = k[i].scale;
        key->theta = k[i].orientation;
        /*  inferring the discrete coordinates in the scale-space grid */
        // We look for the pair of integer (o,s) such that nspo*o+s is the nearest
        // to alpha = nspo * log( k[i].scale / sigma_min) /M_LN2
        // with the constraint s=1,2,..,nspo.
        int o,s;
        int a = (int)(round( n_spo * log( k[i].scale / sigma_min) /M_LN2  ));
        o = (a-1)/n_spo;
        if (o > -1){
            s = (a-1)%n_spo + 1;
        }
        else{
            o = 0;
            s = 0;
        }
        key->o = o;
        key->s = s;
        key->i = (int)( key->x / ( delta_min * exp( key->o * M_LN2)) + 0.5 );
        key->j = (int)( key->y / ( delta_min * exp( key->o * M_LN2)) + 0.5 );
        sift_add_keypoint_to_list(key,keys);
    }
    return keys;
}












/** @brief Main SIFT routine
 * 
 * takes one image as input.        
 * outputs the extracted keypoints in the standard output.
 * 
 * @param flag for the SIFT transform
 * 
 * 0 = one image -> one txt file 
 * 1 = one image -> all txt files
 * 2 = one image -> all txt files + scalespace and DoG
 * 
 */
int main(int argc, char **argv)
{
    // Setting default parameters
    struct sift_parameters* p = sift_assign_default_parameters();
    int flagverb_keys = 0;
    int flagverb_ss = 0;
    char label_ss[256];
    char label_keys[256];
    strcpy(label_ss, "extra");
    strcpy(label_keys, "extra");
    // EXTRA DENSE
    int flag_semigroup = 1;
    int flag_dct = 0;
    int flag_log = 0;
    int flag_interp = 0;

    // Parsing command line
    //int res = parse_options(argc, argv, p, &flagverb_keys, &flagverb_ss, label_keys, label_ss);
    int res = parse_options(argc, argv, p, &flagverb_keys, &flagverb_ss, label_keys, label_ss,
                                      &flag_semigroup, &flag_dct, &flag_log, &flag_interp);
    if (res == EXIT_FAILURE)
        return EXIT_FAILURE;

    /** Loading image */
    int w, h;
    _myfloat* x = read_image(argv[1], &w, &h);
    if(!x)
        fatal_error("File \"%s\" not found.", argv[1]);


    //////////////////////////////
    /** Read keypoints list */
    int n;
    // PROBLEME LA     ne tient pas compte des parametres de discretisation.
    struct sift_keypoint_std * kstd = sift_read_keyslocation_from_file(argv[2], &n);
    struct sift_keypoints* kana = sift_translate_standard_into_anatomy_with_PARAM(kstd, n, p);




    /** Memory dynamic allocation */
    // WARNING 6 steps of the algorithm are recorded.
    struct sift_keypoints **kk = (sift_keypoints **)xmalloc(6*sizeof(struct sift_keypoints*));
    for(int i = 0; i < 6; i++)
        kk[i] = sift_malloc_keypoints();
    // WARNING 4 scale-space representation are recorded (Gaussian, Laplacian, two gradient components)
    struct sift_scalespace **ss = (sift_scalespace **)xmalloc(4*sizeof(struct sift_scalespace*));


    // HERE PATCH
    _myfloat* patches = sift_anatomy_dense_patch(kana, x, w, h, p, ss, kk, flag_semigroup, flag_dct, flag_log, flag_interp);

    /** OUTPUT */
   // int flag = flagverb_keys + 1;
    print_patches(patches, kana->size);

    char name[FILENAME_MAX];
    if(flagverb_keys == 1){
        sprintf(name,"extra_NES_%s.txt",label_keys);              sift_save_keypoints(kk[0], name, 0);
        sprintf(name,"extra_DoGSoftThresh_%s.txt",label_keys);    sift_save_keypoints(kk[1], name, 0);
        sprintf(name,"extra_ExtrInterp_%s.txt",label_keys);       sift_save_keypoints(kk[2], name, 0);
        sprintf(name,"extra_DoGThresh_%s.txt",label_keys);        sift_save_keypoints(kk[3], name, 0);
        sprintf(name,"extra_OnEdgeResp_%s.txt",label_keys);       sift_save_keypoints(kk[4], name, 0);
        sprintf(name,"extra_FarFromBorder_%s.txt",label_keys);    sift_save_keypoints(kk[5], name, 0);
    }
    if (flagverb_ss == 1){
        sprintf(name,"scalespace_%s",label_ss);     print_sift_scalespace_gray_nearestneighbor(ss[0],name);
        sprintf(name,"DoG_%s",label_ss);            print_sift_scalespace_rgb(ss[1],name);
    }

    /* memory deallocation */
    xfree(x);
    xfree(p);
   // sift_free_keypoints(k);
    for(int i = 0; i < 6; i++){
      //  sift_free_keypoints(kk[i]);   // TODO pour gradual
    }
    xfree(kk);
    for(int i = 0; i < 4; i++){
       // sift_free_scalespace(ss[i]);   // TODO pour gradual
    }
    xfree(ss);
    return EXIT_SUCCESS;
}

