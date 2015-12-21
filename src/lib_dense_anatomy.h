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





void compute_and_write_dct_scalespace(FILE* fp,
                                      _myfloat* in, int w_in, int h_in,
                                      struct sift_parameters* p,
                                      int flag_interp,
                                      int flag_type);

// ext is the value added to the cube in case it's out of the scalespace limit
void crop_from_scalespace_file(_myfloat* cube, FILE* fp,
                               int o, int s, int i, int j,
                               int ri, int rj, int rs, _myfloat ext);






#endif //_LIB_DENSE_ANATOMY_H_
