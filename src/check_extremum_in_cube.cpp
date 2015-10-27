/**
 *
 *  To check that there's an extrema in the middle of the cube
 *
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

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
    fprintf(stderr, "Usage:  check_extremum_in_cube cube_file [options]    \n");
    fprintf(stderr, " -r    , the cube is (2r+1) wide    \n");
    fprintf(stderr, " -epsilon                         \n");
    fprintf(stderr, " -r1                              \n");
    fprintf(stderr, " -r2                              \n");
    fprintf(stderr, " -r3                              \n");
    fprintf(stderr, " -type                            \n");

}


// ICI xmalloc et xfree




/**  @brief Process one element of the command line 
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

/**  @brief Process the command line
 *
 */
static int parse_options(int argc, char** argv,
                         int* r,
                         _myfloat* epsilon,
                         _myfloat* r1,
                         _myfloat* r2,
                         _myfloat* r3,
                         int* type)
{
    int isfound;
    char val[128];

    // the cube's width is (2r+1)
    isfound = pick_option(&argc, &argv, "r", val);
    if (isfound ==  1)    *r = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    // epsilon, the tolerance for comparing values
    isfound = pick_option(&argc, &argv, "epsilon", val);
    if (isfound ==  1)    *r = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    // criterion definition - r1
    isfound = pick_option(&argc, &argv, "r1", val);
    if (isfound ==  1)    *r1 = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    // criterion definition - r2
    isfound = pick_option(&argc, &argv, "r2", val);
    if (isfound ==  1)    *r2 = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    // criterion definition - r3
    isfound = pick_option(&argc, &argv, "r3", val);
    if (isfound ==  1)    *r3 = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    // criterion type
    isfound = pick_option(&argc, &argv, "type", val);
    if (isfound ==  1)    *type = atof(val);
    if (isfound == -1)    return EXIT_FAILURE;

    // check for unknown option call
    for(int i = 0; i < argc; i++){
        if (argv[i][0] == '-'){
            fprintf(stderr, "Fatal error: option \"-%s\" is unknown.\n", argv[i]+1);
            print_usage();
            return EXIT_FAILURE;
        }
    }
    // check a name is provided for the cube binary file
    if (argc != 2){
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

    char name_cube[FILENAME_MAX];
    strcpy(name_cube, "extra");
    int r = 1;
    _myfloat epsilon = FLT_EPSILON;
    FILE* fp;

    // criteria parameter (balls)
    _myfloat r1 = 0;
    _myfloat r2 = 0.1;
    _myfloat r3 = 1.0;
    int type = 1;

    // Parsing command line
    int res = parse_options(argc, argv, &r, &epsilon, &r1, &r2, &r3, &type);
    if (res == EXIT_FAILURE)
        return EXIT_FAILURE;

    // Load the cube
    int h = 2*r+1;
    _myfloat* cube = (_myfloat*)xmalloc(h*h*h*sizeof(_myfloat));
    fp = fopen(argv[1], "rb");
    if(!fp)
        fatal_error("Cube file \"%s\" not found", argv[1]);
    //    
    fread(cube, sizeof(_myfloat), h*h*h, fp);
    fclose(fp);

    // Check if there's an extremum inside
    _myfloat output;
    output = check_discrete_extremum(cube, r, epsilon);

    // another criteria TEMP
    //debug(" criteria %i %f %f %f ", r, r1, r2, r3);
    output = confirm_extremum_is_present_inside_ball(cube, r, epsilon, r1, r2, r3, type);


    // Output
    fprintf(stdout, "%f\n", output);

    return EXIT_SUCCESS;
}





