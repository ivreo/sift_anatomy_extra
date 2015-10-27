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
    fprintf(stderr, "Usage:  gradual_sift_extrema image [options...]                                        \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -ss_noct        (8)  number of octaves                                  \n");
    fprintf(stderr, "   -ss_nspo        (3)  number of scales per octaves                       \n");
    fprintf(stderr, "   -ss_dmin      (0.5)  the sampling distance in the first octave          \n");
    fprintf(stderr, "   -ss_smin      (0.8)  blur level on the seed image                       \n");
    fprintf(stderr, "   -ss_sin       (0.5)  assumed level of blur in the input image           \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -thresh_dog (0.0133) threshold over the DoG response                    \n");
    fprintf(stderr, "   -thresh_edge   (10)  threshold over the ratio of principal curvature    \n");
    fprintf(stderr, "         ########### NO DESCRIPTION  ###################                   \n");
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
    fprintf(stderr, "   -flag_semigroup   BOOL (0)    semigroup (1) or direct (0)               \n");
    fprintf(stderr, "   -flag_dct         BOOL (1)    dct (1) or discrete (0)                   \n");
    fprintf(stderr, "   -flag_log         BOOL (0)    normalized Laplacian (1) or DoG (0)       \n");
    fprintf(stderr, "   -flag_interp     (3)  bilin (0) / DCT (1)/ bsplines (3,5,..,11)         \n");
    fprintf(stderr, "   -itermax             5        max number of iterations                  \n");
    fprintf(stderr, "          NOTE: if itermax=0 then all the discrete extrema are output      \n");
    fprintf(stderr, "                              no threshold is applied                      \n");
    fprintf(stderr, "   -epsilon        FLT_EPSILON    (for _myfloat comparison)                \n");
    fprintf(stderr, "   -dog_nspo (int, 3) number of scales per octaves for DoG definition      \n");
    fprintf(stderr, "   -ofstMax_X   (0.5)   interpolation validity domain definition in space  \n");
    fprintf(stderr, "   -ofstMax_S   (0.5)                 ... in scale                         \n");
    fprintf(stderr, "   -flag_jumpinscale   (0 in gradual / not an option in gradual)           \n");
    fprintf(stderr, "   -discrete_extrema_dR (1) The sample is compared to the samples on the surface  \n");
    fprintf(stderr, "                      of a volume containing (2 x halfR + 1)^3 samples     \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -discrete_sphere_r (1) The sample is compared to the samples on         \n");
    fprintf(stderr, "                          the surface of a sphere  x halfR + 1)^3 samples  \n");
    fprintf(stderr, "   -discrete_sphere_dr (1)                                                 \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -extrema_file  name   filename   the keypoints to be tested             \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -outTested  file where the extrema actually tested are stored (border effect) \n");
    fprintf(stderr, "   -outExtrema          filename                                           \n");
    fprintf(stderr, "   -outNotExtrema       filename                                           \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -ball_r1      |                                                         \n");
    fprintf(stderr, "   -ball_r2      |-- define the middle and the surface of the ball to      \n");
    fprintf(stderr, "   -ball_r3      |             confirm the presence of an extremum         \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "   -ball_alpha  |__ a normalized version  (NOTE: not used)                 \n");
    fprintf(stderr, "   -ball_beta   |                                                          \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "                                                                           \n");

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
                         int *flag_interp,
                         char* extrema_file,
                         char* outTested,
                         char* outExtrema,
                         char* outNotExtrema) //  filename of the list of keypoints to test.
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
        p->n_spo = atoi(val);   // redundant
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

    // controlling the extraction of 3d discrete extrema
    isfound = pick_option(&argc, &argv, "discrete_extrema_dR", val);
    if (isfound ==  1)    p->discrete_extrema_dR = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "discrete_sphere_r", val);
    if (isfound ==  1)    p->discrete_sphere_r = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "discrete_sphere_dr", val);
    if (isfound ==  1)    p->discrete_sphere_dr = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;
    
    // the filename of input keypoint file
    isfound = pick_option(&argc, &argv, "extrema_file", val);
    if (isfound ==  1){
        strcpy(extrema_file, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;
    // the filename where the extrema actually tested are stored
    isfound = pick_option(&argc, &argv, "outTested", val);
    if (isfound ==  1){
        strcpy(outTested, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;
    // the filename where the extrema that are confirmed to be extrema
    isfound = pick_option(&argc, &argv, "outExtrema", val);
    if (isfound ==  1){
        strcpy(outExtrema, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;
    // the filename where the extrema that are confirmed to be extrema
    isfound = pick_option(&argc, &argv, "outNotExtrema", val);
    if (isfound ==  1){
        strcpy(outNotExtrema, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;

    // NOT USED TODO REMOVE
    isfound = pick_option(&argc, &argv, "ball_r1", val);
    if (isfound ==  1)    p->ball_r1 = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "ball_r2", val);
    if (isfound ==  1)    p->ball_r2 = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "ball_r3", val);
    if (isfound ==  1)    p->ball_r3 = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;
    // NORMALIZED BALL CRITERION
    isfound = pick_option(&argc, &argv, "ball_alpha", val);
    if (isfound ==  1)    p->ball_alpha = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;
    isfound = pick_option(&argc, &argv, "ball_beta", val);
    if (isfound ==  1)    p->ball_beta = atof(val);
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
   // float k_nspo =  exp( M_LN2/( (float)(dog_nspo)));
   // float k_3 =  exp( M_LN2/( (float)3));
    float k_nspo = pow(2, (float)(dog_nspo));
    float k_3 =  pow(2, 3.0);
    float factor = (k_3 - 1) / (k_nspo - 1);

    for (int i = 0; i < keys->size; i++){
        struct keypoint* k = keys->list[i];
        float val = k->val * factor;
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
        // k->s is always equal to 1 for the 'gradual' implementation // TODO correct the structure somewhere before output.
        fprintf(stdout, "%i %i %i %f %f %f ", k->i, k->j, k->s,  k->ofstX,  k->ofstY,  k->ofstS);
        fprintf(stdout, "\n");
    }
}


// TODO - factorize code with the print_keypoints_and_vals function
static void save_keypoints_and_vals(const char* filename, const struct sift_keypoints* keys, int dog_nspo)
{

    FILE* file = fopen(filename, "w");

   // float k_nspo =  exp( M_LN2/( (float)(dog_nspo)));
   // float k_3 =  exp( M_LN2/( (float)3));
    float k_nspo = pow(2, (float)(dog_nspo));
    float k_3 =  pow(2, 3.0);
    float factor = (k_3 - 1) / (k_nspo - 1);

    for (int i = 0; i < keys->size; i++){
        struct keypoint* k = keys->list[i];
        float val = k->val * factor;
#ifdef QUAD // conversion to string and then fprintf
        char str[1024]; // long strinpl
        quadmath_snprintf(str, sizeof(str), "%.30Qg %.30Qg %.30Qg %.30Qg %.30Qg",  k->x, k->y, k->sigma, val);
        fprintf(file, "%s \n", str);
#else
        fprintf(file, "%f %f %f %f ", k->x, k->y, k->sigma, val);
#endif
        // The 27 points in the cube.
        for(int n = 0; n < 27; n++){
            fprintf(file, "%33.30f ", k->neighbors[n]);
        }
        // EXTRA - the interpolation offset for the last interpolation and grid position
        // k->s is always equal to 1 for the 'gradual' implementation // TODO correct the structure somewhere before output.
        fprintf(file, "%i %i %i %f %f %f ", k->i, k->j, k->s,  k->ofstX,  k->ofstY,  k->ofstS);
        fprintf(file, "\n");
    }


    fclose(file);
}






//TODO -- TO BE TESTED THOROUGHLY -- DOUBLE CHECK
static void find_nearest_scalespace_sample(_myfloat x,     // keypoint position
                                    _myfloat y,
                                    _myfloat sigma,
                                    int n_spo,             // scalespace parameters
                                    _myfloat sigma_min,
                                    _myfloat delta_min,
                                    int *o,          // nearest sample coordinate
                                    int *s,
                                    int *i,
                                    int *j)
{
    int a = (int)(round( n_spo * log( sigma / sigma_min) /M_LN2  ));
    *o = (a-1)/n_spo;
    if (*o > -1){
        *s = (a-1)%n_spo + 1;
    }
    else{
        *o = 0;
        *s = 0;
    }
    *i = (int)( x / ( delta_min * exp( *o * M_LN2)) + 0.5 );
    *j = (int)( y / ( delta_min * exp( *o * M_LN2)) + 0.5 );
}
//TODO: A simple way to visualize if this works is to read interpolated
//keypoints with a very dense scalespace, save them as discrete keypoints and
//then map the interpolated keypoints against the discrete keypoints



void read_keypoints_and_vals(struct sift_keypoints* keys,
                             const char* name,
                             int n_spo,
                             _myfloat sigma_min,
                             _myfloat delta_min)
{
    size_t buffer_size = 1024 * 1024;  // 1MB buffer for long lines.
    char* buffer = (char*)xmalloc(buffer_size);   // note: we typecast malloc (this is C++).
    FILE* stream = fopen(name,"r");
    if ( !stream)
        fatal_error("File \"%s\" not found.", name);
    while(fgets(buffer, buffer_size, stream) != NULL){
        int pos = 0;
        int read = 0;
        struct keypoint* key = sift_malloc_keypoint(8, 4, 36); // n_ori, n_hist, n_bins);

        // read coordinates
        float fx, fy, fsigma;
        sscanf(buffer+pos,"%f  %f  %f  %n", &fx
                                          , &fy
                                          , &fsigma
                                          , &read);
        pos+=read;
   //     fprintf(stderr, " reading %f %f %f  \n", fx, fy, fsigma); 

        // datatype convertion
        _myfloat x = (_myfloat)fx;
        _myfloat y = (_myfloat)fy;
        _myfloat sigma = (_myfloat)fsigma;
        
        // find the nearest scalespace sample
        int o, s, i ,j;
        find_nearest_scalespace_sample(x,y,sigma,                   // keypoint coordinates
                                       n_spo, sigma_min, delta_min, // scalespace parameters  
                                       &o, &s, &i, &j);             // (output) coordinates of nearest sample.

        // saving keypoint
        key->x = x;
        key->y = y;
        key->sigma = sigma;
        key->o = o;
        key->s = s;
        key->i = i;
        key->j = j;

        // DEBUG
        //fprintf(stderr, "DEBUGin_read_keypoints_and_vals: Candidate keypoint (s,i,j) = (%i,%i,%i) \n", s,i,j);

        sift_add_keypoint_to_list(key, keys);

        // read the rest of the line
        // (to go to the next line)
        // (MAYBE there's a clever way)
        //
        // TODO - read the realdata instead - can be useful to keep the information of the surrounding !!!
        //
        for(int i = 0; i < (27+6); i++){
            float tmp;
            sscanf(buffer+pos, "%f %n", &tmp, &read);
            pos +=read;
        }
    }
    xfree(buffer);
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
    int flag_semigroup = 0; // This has no effect anyway - everything is computed from the input image after interpolation
    int flag_dct = 1;
    int flag_log = 0;
    int flag_interp = 3;
    // name of input file
    char extrema_file[256];   // TODO added proper test to handle 'file doesn't exist' situation
    char outTested[256];      // TODO added proper test to handle 'file doesn't exist' situation
    char outExtrema[256];     // TODO added proper test to handle 'file doesn't exist' situation
    char outNotExtrema[256];  // TODO added proper test to handle 'file doesn't exist' situation
    
    // Parsing command line
    int res = parse_options(argc, argv, p, &flagverb_keys, &flagverb_ss, label_keys, label_ss,
                                      &flag_semigroup, &flag_dct, &flag_log, &flag_interp,
                                      extrema_file, outTested, outExtrema, outNotExtrema);// check_extrema
    if (res == EXIT_FAILURE)
        return EXIT_FAILURE;

    // Loading image
    int w, h;
    _myfloat* x = read_image(argv[1], &w, &h);
    if(!x)
        fatal_error("File \"%s\" not found.", argv[1]);


    // Loading list of keypoints
    struct sift_keypoints* keysIn = sift_malloc_keypoints();
    //fprintf(stderr, "     number of scales considered in the conversion to discrete sample %i \n",  p->n_spo);
    read_keypoints_and_vals(keysIn, extrema_file, p->n_spo, p->sigma_min, p->delta_min);

    // Checking for extrema
    struct sift_keypoints* keysTested = sift_malloc_keypoints();
    struct sift_keypoints* keysExtrema = sift_malloc_keypoints();
    struct sift_keypoints* keysNotExtrema = sift_malloc_keypoints();
    sift_anatomy_gradual_check_extrema(x, w, h,
                                       keysIn,
                                       keysTested,
                                       keysExtrema,
                                       keysNotExtrema,
                                       p,
                                       flag_semigroup,
                                       flag_dct,
                                       flag_interp);


    // TMP
    fprintf(stderr, "size input %i  \n", keysIn->size);
    fprintf(stderr, "size tested %i  \n", keysTested->size);
    fprintf(stderr, "size Extrema %i  \n", keysExtrema->size);
    fprintf(stderr, "size Not Extrema %i  \n", keysNotExtrema->size);

    //print_keypoints_and_vals(keysExtrema, p->dog_nspo);
    save_keypoints_and_vals(outTested, keysTested, p->dog_nspo);
    save_keypoints_and_vals(outExtrema, keysExtrema, p->dog_nspo);
    save_keypoints_and_vals(outNotExtrema, keysNotExtrema, p->dog_nspo);

    xfree(x);
    xfree(p);
    sift_free_keypoints(keysTested);
    sift_free_keypoints(keysExtrema);
    sift_free_keypoints(keysNotExtrema);

    return EXIT_SUCCESS;
}




