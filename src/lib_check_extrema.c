/** @file lib_check_extrema.c
 *  20151023
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
//#include <math.h>

#include "lib_scalespace.h"
#include "lib_keypoint.h"
#include "lib_util.h"


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



void read_keypoints(struct sift_keypoints* keys,
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
        struct keypoint* key = sift_malloc_keypoint(8, 4, 36); // n_ori, n_hist, n_bins);

        // read coordinates
        float fx, fy, fsigma;
        sscanf(buffer,"%f  %f  %f%*[^\n]", &fx, &fy, &fsigma);

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
        sift_add_keypoint_to_list(key, keys);
    }
    xfree(buffer);
}


/** @brief Extract from the scale-space
 *         a cube of samples around the discrete
 *         sample
 *
 *
 */
void extract_portion_of_scalespace(_myfloat* cube, struct sift_scalespace*d,  int r, int o, int s, int i, int j)
{
    _myfloat* imStack = d->octaves[o]->imStack;
    int h = d->octaves[o]->h;
    int w = d->octaves[o]->w;
    int nSca = d->octaves[o]->nSca;
    int n = 0;
    for(int ks = s-r; ks <= s+r; ks++){
        for(int ki = i-r; ki <= i+r; ki++){
            for(int kj = j-r; kj <= j+r; kj++){
                if ( (ks>=0) && (ki>=0) && (kj>=0) && (ks<nSca) && (ki<h) && (kj<w) ){
                    cube[n] = imStack[ks*w*h + ki*w + kj];
                }else{
                    cube[n] = NAN;
                }
                n++;
            }
        }
    }
}



/** @brief  Compare a sample to its 26 neighbors 
 *
 *  Output:
 *    1   - discrete maximum
 *   -1   - discrete minimum
 *    0   - nothing
 *    NAN - we went beyond the limits of the scale-space 
 */
_myfloat check_discrete_extremum(_myfloat* cube, int r, _myfloat epsilon)
{
    _myfloat output = 0.;

    int h = (2*r+1);
    int middle = (h*h*h-1)/2;
    _myfloat center_value = cube[middle];

    // Positions of neighbors
    int neighbors[h*h*h];
    int n = 0;
    for (int ds = -1; ds <= 1; ds++) {
        for (int di = -1; di <= 1; di++) {
            for (int dj = -1; dj <= 1; dj++) {
                if (ds != 0 || di != 0 || dj != 0) {
                    neighbors[n] = middle + (ds * h + di) * h + dj;
                    n++;
                }
            }
        }
    }
    // Gather neighbors values
    _myfloat values[h*h*h];
    bool is_there_nan = false;
    for (int l = 0; l < n; l++){
        values[l] = cube[neighbors[l]];
        if (isnan(values[l])){
            is_there_nan = true;
        }
    }
    // THE TEST
    // Checking the neighbors are not outside the limits of the scale-space
    if (is_there_nan == true){
        output = NAN;
    } else {
        bool is_local_min = true;
        // An optimizing compiler will unroll this loop.
        for (int l = 0; l < n; l++) {
            if (values[l] - epsilon <= center_value) {
                is_local_min = false;
                break; // Can stop early if a smaller neighbor was found.
            }
        }
        bool is_local_max = true;
        // Can skip max check if center point was determined to be a local min.
        if (is_local_min) {
            is_local_max = false;
        } else {
            // An optimizing compiler will unroll this loop.
            for (int l = 0; l < n; l++) {
                if (values[l] + epsilon >= center_value) {
                    is_local_max = false;
                    break; // Can stop early if a larger neighbor was found.
                }
            }
        }
        if (is_local_max == true)
            output = 1.;
        if (is_local_min == true)
            output = -1.;
    }
    // END - THE TEST
    return output;
}




/**
 *  @brief  Check there's a extremum in the middle of the cube
 *
 *
 *  Output:
 *    1   - discrete maximum
 *   -1   - discrete minimum
 *    0   - nothing
 *    NAN - we went beyond the limits of the scale-space 
 *
 *  Type:
 *    0 - classic criteria
 *    1 - mean value on the center
 *    2 - JMM (max in > max_out) !!
 *
 */
_myfloat confirm_extremum_is_present_inside_ball(_myfloat* cube,
                                                 int r,
                                                 _myfloat epsilon,
                                                 _myfloat r1,
                                                 _myfloat r2,
                                                 _myfloat r3,
                                                 int type)
{
    _myfloat output = 0.;
    int h = (2*r+1);
    int middle = (h*h*h-1)/2;

    // Checking that the cube is large enough
    if( ceil(r3)>r)
        fatal_error(" cube too little (r=%i) for the criteria (r3=%f)", r, r3);

    // Positions of the neighbors in the grid
    int neighbor_in[h*h*h];
    int neighbor_out[h*h*h];
    int n_in = 0;
    int n_out = 0;
    for (int ds = -r; ds <= r; ds++) {
        for (int di = -r; di <= r; di++) {
            for (int dj = -r; dj <= r; dj++) {
                // compute the sample distance to the candidate extrema
                _myfloat dist2 = di*di + dj*dj + ds*ds;
                // list of neighbors that are in the middle of the  ball
                if (dist2 < r1*r1) {
                    neighbor_in[n_in] = middle + (ds * h + di) * h + dj;
                    n_in++;
                }
                // list of neighbors that are on the surface of the ball
                if ((dist2 >= r2*r2) && (dist2 < r3*r3)) {
                    neighbor_out[n_out] = middle + (ds * h + di) * h + dj;
                    n_out++;
                }
            }
        }
    }
    debug( "n_in=%i  n_out=%i", n_in, n_out);

    // Load values
    _myfloat values_in[h*h*h];
    _myfloat values_out[h*h*h];
    // - in
    for (int n = 0; n < n_in; n++){
        values_in[n] = cube[neighbor_in[n]];
    }
    // - out
    bool is_there_nan = false;
    for (int n = 0; n < n_out; n++){
        values_out[n] = cube[neighbor_out[n]];
        if( isnan(values_out[n])){
            is_there_nan = true;
        }
    }
    // Compute mins, maxes, means and stds
    _myfloat max_in =  array_max(values_in, n_in);
    _myfloat min_in =  array_min(values_in, n_in);
    _myfloat max_out =  array_max(values_out, n_out);
    _myfloat min_out =  array_min(values_out, n_out);
    _myfloat mean_in, std_in;
    _myfloat mean_out, std_out;
    array_mean_and_std(values_in, n_in, &mean_in, &std_in);
    array_mean_and_std(values_out, n_out, &mean_out, &std_out);

    // THE TEST
    // Checking the neighbors are not outside the limits of the scale-space
    if (is_there_nan == true){
        output = NAN;
    } else {
        bool is_local_max;
        bool is_local_min;
        if (type == 0){
            // 
            is_local_max = (min_in > max_out + epsilon);
            is_local_min = (max_in < min_out - epsilon);
        }else if(type == 1){
            // Denoised criterion
            // (this is not exactly Mauricio's since the values on the surface are not smoothed *locally*).
            is_local_max = (mean_in > max_out + epsilon);
            is_local_min = (mean_in < min_out - epsilon);
        }else{
            // JM criterion
            is_local_max = (max_in > max_out);
            is_local_min = (min_in < min_out);
        }


        if (is_local_max == true)
            output = 1.;
        if (is_local_min == true)
            output = -1.;
    }
    // END - THE TEST
    return output;
}







