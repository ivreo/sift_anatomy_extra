/*
IPOL SIFT
Copyright (C) 2014, Ives Rey Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>
 */



#ifndef _LIB_DENSE_ANATOMY_H_
#define _LIB_DENSE_ANATOMY_H_


struct sift_keypoints* sift_anatomy_dense(_myfloat* x, int w, int h, struct sift_parameters* p,
                                          struct sift_scalespace* ss[4],
                                          struct sift_keypoints* kk[6],
                                          int flag_interp);

struct sift_keypoints* sift_anatomy_gradual_old(_myfloat* x, int w, int h, struct sift_parameters* p,
                                          struct sift_keypoints* kk[6],
                                          int flag_dct,
                                          int flag_interp);

struct sift_keypoints* sift_anatomy_gradual(_myfloat* x, int w, int h, struct sift_parameters* p,
                                          struct sift_keypoints* kk[6],
                                          int flag_interp);




#endif //_LIB_DENSE_ANATOMY_H_
