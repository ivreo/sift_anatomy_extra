#ifndef _LIB_ELLIPSE_H_
#define _LIB_ELLIPSE_H_

#define N_BINS 36 // for orientation attribution
#define N_ORI 8
#define N_HIST 4
#define N_DESCR 128
#define DESCR_THRESH 0.2



#include "lib_keypoint.h"


#include "lib_sift_anatomy.h" // FIXME  TODO for the definition of the parameter // for match_cli_ellipse only (that doesn't include lib_sift_anatomy.h)

void sift_read_ellipses(struct sift_keypoints* keys,
                        const char* name,
                        int n_hist,
                        int n_ori,
                        int n_bins,
                        int flag);

void fprintf_one_ellipse(FILE* f, const struct keypoint* k, int n_descr, int n_bins, int flag);

void fprintf_ellipses(FILE* f, const struct sift_keypoints* keys, int flag);

void sift_save_ellipses(const struct sift_keypoints* keys, const char* name, int flag);

void sift_print_ellipses(const struct sift_keypoints* keys, int flag);

void print_pairs_ellipses(const struct sift_keypoints *k1,
                          const struct sift_keypoints *k2);

void ellipse_sift_anatomy_orientation_and_description(const _myfloat* x,
                                              int w,
                                              int h,
                                              const struct sift_parameters* p,
                                              struct sift_keypoints* k, int gauss_flag);

void estimate_scalespace_index(struct sift_keypoints *l, struct sift_parameters* p);

#endif /* _LIB_ELLIPSE_H_ */
