/**
 *  To extract neighbors from a keypoint
 *
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

extern "C" {

#include "lib_sift_anatomy.h"
#include "lib_check_extrema.h"
//#include "lib_io_scalespace.h"
//#include "io_png.h"
//#include "lib_util.h"
//#include "iio.h"
//#include "lib_dense_anatomy.h"

}
//#include "io_exr.h"




void print_usage()
{
    fprintf(stderr, "Usage:  grab_ss_neighborhood \n");
    fprintf(stderr, " -dog   [name scale-space]                                                 \n");
    fprintf(stderr, " -keys [name keypoints]                                                    \n");
    fprintf(stderr, " -r                    (the cube is (2r+1) wide                            \n");
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
                         char* label_keys,
                         char* label_ss,
                         int* r)
{
    int isfound;
    char val[128];


    // keypoints file
    isfound = pick_option(&argc, &argv, "keys", val);
    if (isfound ==  1){
        strcpy(label_keys, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;

    // scalespace file
    isfound = pick_option(&argc, &argv, "dog", val);
    if (isfound ==  1){
        strcpy(label_ss, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;

    // dimension of the cube
    isfound = pick_option(&argc, &argv, "r", val);
    if (isfound ==  1)    *r = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    // check for unknown option call
    for(int i = 0; i < argc; i++){
        if (argv[i][0] == '-'){
            fprintf(stderr, "Fatal error: option \"-%s\" is unknown.\n", argv[i]+1);
            print_usage();
            return EXIT_FAILURE;
        }
    }
    //
    if (argc != 1){
        print_usage();
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}









/** @brief Main
 *
 *
 */
int main(int argc, char **argv)
{
    char name_keys[FILENAME_MAX];
    char name_dog[FILENAME_MAX];
    char name[FILENAME_MAX];
    strcpy(name_dog, "extra");
    strcpy(name_keys, "extra");
    strcpy(name, "extra");
    int r = 2;
    FILE* fp;
    
    // Parsing command line
    int res = parse_options(argc, argv, name_keys, name_dog, &r);
    if (res == EXIT_FAILURE)
        return EXIT_FAILURE;

    // Loading the DoG scale-space
    fp = fopen(name_dog, "rb");
    struct sift_scalespace* d = sift_read_scalespace_binary_file(fp);
    fclose(fp);

    // guessing parameters from the loaded DoG scale-space
    int n_spo = d->octaves[0]->nSca - 2; // only 2 extra scales in the image stack
    _myfloat sigma_min = d->octaves[0]->sigmas[0];
    _myfloat delta_min = d->octaves[0]->delta;

    debug("nspo %i", n_spo);

    // Loading list of keypoints
    struct sift_keypoints* keysIn = sift_malloc_keypoints();
    read_keypoints(keysIn, name_keys, n_spo, sigma_min, delta_min);


    // Extracting the cube around each keypoint
    int h = (2*r+1);
    _myfloat* cube = (_myfloat*)xmalloc(h*h*h*sizeof(_myfloat));
    for(int k = 0; k < keysIn->size; k++){

        int o = keysIn->list[k]->o;
        int s = keysIn->list[k]->s;
        int i = keysIn->list[k]->i;
        int j = keysIn->list[k]->j;

        extract_portion_of_scalespace(cube, d, r, o, s, i, j);
        sprintf(name, "%s_local_%05i_outof_%05i.bin", name_keys, k, keysIn->size);
        fp = fopen(name, "wb");
        fwrite(cube, sizeof(_myfloat), h*h*h, fp);
        fclose(fp);
    }
    return EXIT_SUCCESS;
}


