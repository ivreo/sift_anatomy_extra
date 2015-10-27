/** @file lib_check_extrema.c
 *  20151023
 *
 *
 */


#include "lib_util.h" // for _myfloat definition
#include "lib_scalespace.h"


#ifndef _LIB_CHECK_EXTREMA_H_
#define _LIB_CHECK_EXTREMA_H_


/** @brief Read keypoints from a file
 *
 */
void read_keypoints(struct sift_keypoints* keys,
                    const char* name,
                    int n_spo,
                    _myfloat sigma_min,
                    _myfloat delta_min);


/** @brief Extract sample cube  from scale-space
 *
 */
void extract_portion_of_scalespace(_myfloat* cube,
                                   struct sift_scalespace* d,
                                   int r,
                                   int o,
                                   int s,
                                   int i,
                                   int j);


/** @brief  Compare a sample to its 26 neighbors 
 *
 *  Output:
 *    1   - discrete maximum
 *   -1   - discrete minimum
 *    0   - nothing
 *    NAN - we went beyond the limits of the scale-space 
 */
_myfloat check_discrete_extremum(_myfloat* cube, int r,
                                 _myfloat epsilon);


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
                                                 int type);




#endif
