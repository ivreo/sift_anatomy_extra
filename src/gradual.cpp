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

}
#include "io_exr.h"




void print_usage()
{
    fprintf(stderr, "Anatomy of the SIFT method (www.ipol.im/pub/pre/82/)  ver 20140801         \n");
    fprintf(stderr, "Usage:  sift_cli image [options...]                                        \n");
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
    fprintf(stderr, "   -ss_fnspo (3.00) float number of scales per octaves (-ss_nspo not used) \n");
    fprintf(stderr, "   -flag_interp     (0)  bilin (0) / DCT (1)/ bsplines (3,5,..,11)         \n");
    fprintf(stderr, "   -itermax             5        max number of iterations                  \n");
    fprintf(stderr, "          NOTE: if itermax=0 then all the discrete extrema are output      \n");
    fprintf(stderr, "                              no threshold is applied                      \n");
    fprintf(stderr, "   -epsilon        FLT_EPSILON    (for _myfloat comparison)                \n");
    fprintf(stderr, "   -dog_nspo (int, 3) number of scales per octaves for DoG definition      \n");
    fprintf(stderr, "   -ofstMax_X   (0.5)   interpolation validity domain definition in space  \n");
    fprintf(stderr, "   -ofstMax_S   (0.5)                 ... in scale                         \n");
    fprintf(stderr, "   -flag_jumpinscale   (0 in gradual / not an option in gradual)           \n");
}


//    /**  @brief  The DMIN allowing a balanced scale-space corresponding to NSPO, SSMIN
//     *
//     *
//     */
//    static _myfloat ideal_dmin(int nspo, _myfloat ssmin)
//    {
//        _myfloat dmin = 1/sqrt(2)*ssmin*(1 + pow( 2,(1./nspo)) - pow(2, (1-1./nspo)) );
//        return dmin;
//    }



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
                         int *flag_interp)
{
    int isfound;
    char val[128];

    isfound = pick_option(&argc, &argv, "ss_noct", val);
    if (isfound ==  1)    p->n_oct = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ss_nspo", val);
    if (isfound ==  1){
        p->n_spo = atoi(val);
        p->fnspo = atof(val); // redundant
    }
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
//    isfound = pick_option(&argc, &argv, "flag_semigroup", val);
//    if (isfound ==  1)    *flag_semigroup = atoi(val);
//    if (isfound == -1)    return EXIT_FAILURE;
//    isfound = pick_option(&argc, &argv, "flag_dct", val);
//    if (isfound ==  1)    *flag_dct = atoi(val);
//    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "itermax", val);
    if (isfound ==  1)    p->itermax = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "dog_nspo", val);
    if (isfound ==  1)    p->dog_nspo = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "epsilon", val);
    if (isfound ==  1)    p->epsilon = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "flag_interp", val);
    if (isfound ==  1)    *flag_interp = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "ss_fnspo", val);
    if (isfound ==  1){
        p->fnspo = atof(val);
        p->n_spo = atoi(val); // redundant
    }
    if (isfound == -1)    return EXIT_FAILURE;

    //  controlling the interpolation
    isfound = pick_option(&argc, &argv, "flag_jumpinscale", val);
    if (isfound ==  1)    p->flag_jumpinscale = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "ofstMax_X", val);
    if (isfound ==  1)    p->ofstMax_X = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "ofstMax_S", val);
    if (isfound ==  1)    p->ofstMax_S = atof(val);
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
    if (argc != 2){
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
    else if ( string_suffix(n, ".tif")  || string_suffix(n, ".tiff")  || string_suffix(n, ".TIFF")){
        // IIO
        x = iio_read_image_float(n, w, h); 
        for(int i=0; i < (*w)*(*h); i++)
            //x[i] /= 256.0;
            x[i] /= 4096.0;
    }
    else{
        fatal_error("This version only recognizes png and exr(32bits) images");
    }

    _myfloat* y = (_myfloat*)malloc((*w)*(*h)*sizeof(_myfloat));
    for (int i = 0; i < (*w)*(*h); i++){
        y[i] = (_myfloat)x[i];
    }
    free(x);
    return(y);
}



// Print function
// Note that the value is normalized (equivalent to dog_nspo = 3)
// Note that the value in the 27 long vector are not normalized
void print_keypoints_and_vals(const struct sift_keypoints* keys, int dog_nspo)
{
    _myfloat k_nspo = pow(2, (_myfloat)(dog_nspo));
    _myfloat k_3 =  pow(2, 3.0);
    _myfloat factor = (k_3 - 1) / (k_nspo - 1);

    for (int i = 0; i < keys->size; i++){
        struct keypoint* k = keys->list[i];
        _myfloat val = k->val * factor;
#ifdef QUAD // conversion to string and then fprintf
        char str[1024]; // long strinpl
        quadmath_snprintf(str, sizeof(str), "%.30Qg %.30Qg %.30Qg %.30Qg %.30Qg",  k->x, k->y, k->sigma, val);
        fprintf(stdout, "%s \n", str);
#else
        fprintf(stdout, "%f %f %f %f ", k->x, k->y, k->sigma, val);
#endif
        // The 27 points in the cube.
        for(int n = 0; n < 27; n++){
            fprintf(stdout, "%33.30f ", k->neighbors[n]);
        }
        // EXTRA - the interpolation offset for the last interpolation and grid position
        fprintf(stdout, "%i %i %i %f %f %f ", k->i, k->j, k->s,  k->ofstX,  k->ofstY,  k->ofstS);
        //
        fprintf(stdout, "\n");
    }
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
    int flag_interp = 0;
    // Parsing command line
    int res = parse_options(argc, argv, p, &flagverb_keys, &flagverb_ss, label_keys, label_ss, &flag_interp);
    if (res == EXIT_FAILURE)
        return EXIT_FAILURE;

    /** Loading image */
    int w, h;
    _myfloat* x = read_image(argv[1], &w, &h);
    if(!x)
        fatal_error("File \"%s\" not found.", argv[1]);


    /** Memory dynamic allocation */
    // WARNING 6 steps of the algorithm are recorded.
    struct sift_keypoints **kk = (sift_keypoints **)xmalloc(6*sizeof(struct sift_keypoints*));
    for(int i = 0; i < 6; i++)
        kk[i] = sift_malloc_keypoints();

    /** Algorithm */
//    struct sift_keypoints* k = sift_anatomy_gradual_old(x, w, h, p, kk, flag_dct, flag_interp);
    struct sift_keypoints* k = sift_anatomy_gradual(x, w, h, p, kk, flag_interp);

    /** OUTPUT */
    int flag = flagverb_keys + 1;
    sift_print_keypoints(k, flag);
    //print_keypoints_and_vals(k, p->dog_nspo);

    char name[FILENAME_MAX];
    if(flagverb_keys == 1){
        sprintf(name,"extra_NES_%s.txt",label_keys);              sift_save_keypoints(kk[0], name, 0);
        sprintf(name,"extra_DoGSoftThresh_%s.txt",label_keys);    sift_save_keypoints(kk[1], name, 0);
        sprintf(name,"extra_ExtrInterp_%s.txt",label_keys);       sift_save_keypoints(kk[2], name, 0);
        sprintf(name,"extra_DoGThresh_%s.txt",label_keys);        sift_save_keypoints(kk[3], name, 0);
        sprintf(name,"extra_OnEdgeResp_%s.txt",label_keys);       sift_save_keypoints(kk[4], name, 0);
        sprintf(name,"extra_FarFromBorder_%s.txt",label_keys);    sift_save_keypoints(kk[5], name, 0);
    }

    /* memory deallocation */
    xfree(x);
    xfree(p);
    sift_free_keypoints(k);
    for(int i = 0; i < 6; i++){
        sift_free_keypoints(kk[i]);   // TODO pour gradual
    }
    xfree(kk);


    return EXIT_SUCCESS;
}

