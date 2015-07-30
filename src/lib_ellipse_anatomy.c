
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "lib_description.h"
#include "lib_discrete.h"
#include "lib_sift_anatomy.h"
#include "lib_util.h"

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )
#define ABS(x) ((x)<0?-(x):(x))

#define EPSILON 0



static void ellipse_cartesian_to_polar(float a, _myfloat b, _myfloat c,
                                _myfloat* alpha, _myfloat* r1, _myfloat* r2)
{
    *alpha = 0;
    if ( ABS(a-c) > 0.000000000000000000001)
        *alpha = atan(b/(a-c))/2;
    _myfloat ct = cos(*alpha);
    _myfloat st = sin(*alpha);
    *r1 = sqrt(1./(a*ct*ct + b*ct*st + c*st*st));
    *r2 = sqrt(1./(a*st*st - b*ct*st + c*ct*ct));
}




static void ellipse_sift_accumulate_orientation_histogram(float x_key, _myfloat y_key,
                                           _myfloat a, _myfloat b, _myfloat c,
                                           const _myfloat* imX, const _myfloat* imY,
                                           int w, int h,
                                           int nbins, _myfloat lambda_ori, _myfloat* hist)
{
    /// Ellipse cartesian to polar convertion
    _myfloat r1, r2, alpha;
    ellipse_cartesian_to_polar(a, b, c, &alpha, &r1, &r2);
    _myfloat sigma_key = sqrt(r1*r2);

    /// Initialize output vector
    for(int i= 0;i<nbins;i++){hist[i] = 0.0;}

    /// Contributing pixels are inside a patch [siMin;siMax] X [sjMin;sjMax] of w 6*lambda_ori*sigma_key (=9*sigma_key) 
    _myfloat R = 3*lambda_ori*sigma_key;
    _myfloat Rp = 3*lambda_ori*M_SQRT2*MAX(r1,r2);
    int siMin = MAX(0, (int)(x_key-Rp+0.5));
    int sjMin = MAX(0, (int)(y_key-Rp+0.5));
    int siMax = MIN((int)(x_key+Rp+0.5), h-1);
    int sjMax = MIN((int)(y_key+Rp+0.5), w-1);


    /// For each pixel inside the patch.
    for(int si = siMin; si <= siMax; si++){
        for(int sj = sjMin; sj <= sjMax; sj++){

            _myfloat X, Y;    // coordinates in the normalized referential
            _myfloat dX, dY;  // gradient in the normalized referential
            X = si - x_key;
            Y = sj - y_key;
            apply_rotation(X, Y, &X, &Y, -alpha);
            X  = X*sigma_key/r1;
            Y  = Y*sigma_key/r2;
            apply_rotation(X, Y, &X, &Y, alpha);

            if (MAX(ABS(X),ABS(Y)) < R){

                // gradient in the normalized referential
                dX = imX[si*w+sj];
                dY = imY[si*w+sj];
                apply_rotation(dX, dY, &dX, &dY, -alpha);
                dX  = dX/sigma_key*r1;
                dY  = dY/sigma_key*r2;
                apply_rotation(dX, dY, &dX, &dY, alpha);

                /// the gradient orientation 
                _myfloat ori = modulus(atan2(dY, dX), 2*M_PI);
                // gradient magnitude with Gaussian weighing
                _myfloat t = lambda_ori*sigma_key;
                _myfloat M = hypot(dX, dY) * exp(-(X*X+Y*Y)/(2*t*t));

                /// Determine the bin index in the circular histogram
                //int gamma = (int)(sOri/(2*M_PI)*(float)nbins+0.5)%nbins;
                int gamma = (int)(ori/(2*M_PI)*nbins+0.5)%nbins;

                /// Add the contribution to the orientation histogram
                hist[gamma] += M;
            }
        }
    }
}




static void ellipse_sift_extract_feature_vector(float x_key, _myfloat y_key,
                                         _myfloat a, _myfloat b, _myfloat c,
                                         _myfloat theta_key,
                                         const _myfloat* imX, const _myfloat* imY,
                                         int w, int h,
                                         int Nhist,
                                         int Nbins,
                                         _myfloat lambda_descr,
                                         int gauss_flag,
                                         _myfloat* descr)
{
    // Initialize descr tab
    for(int i = 0; i < Nhist*Nhist*Nbins; i++){descr[i] = 0.0;}

    /// Ellipse cartesian to polar convertion
    _myfloat r1, r2, alpha;
    ellipse_cartesian_to_polar(a, b, c, &alpha, &r1, &r2);
    _myfloat sigma_key = sqrt(r1*r2);

    /// Contributing pixels are inside a patch [siMin;siMax]X[sjMin;sjMax] of 2*lambda_descr*sigma_key*(nhist+1)/nhist
    _myfloat l = (1+1/(float)Nhist)*lambda_descr;
    _myfloat R = l*sigma_key;
    _myfloat Rp = M_SQRT2*l*MAX(r1,r2);
    int siMin = MAX(0, (int)(x_key - Rp + 0.5));
    int sjMin = MAX(0, (int)(y_key - Rp + 0.5));
    int siMax = MIN((int)(x_key + Rp + 0.5), h-1);
    int sjMax = MIN((int)(y_key + Rp + 0.5), w-1);

    // For each pixel inside the patch.
    for(int si = siMin; si < siMax; si++){
        for(int sj = sjMin; sj < sjMax; sj++){

            // Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
            _myfloat X, Y, dX, dY;
            X = si - x_key;
            Y = sj - y_key;
            apply_rotation(X, Y, &X, &Y, -alpha);
            X = X*sigma_key/r1;
            Y = Y*sigma_key/r2;
            apply_rotation(X, Y, &X, &Y, alpha-theta_key);

            // Does this sample fall inside the descriptor area ?
            if (MAX(ABS(X),ABS(Y)) < R) {
                // gradient in the normalized referential
                dX = imX[si*w+sj];
                dY = imY[si*w+sj];
                apply_rotation(dX, dY, &dX, &dY, -alpha);
                dX = dX/sigma_key*r1;
                dY = dY/sigma_key*r2;
                apply_rotation(dX, dY, &dX, &dY, alpha);

                // Compute the gradient orientation (theta) on keypoint referential.
                _myfloat ori = atan2(dY,dX) - theta_key;
                ori = modulus(ori, 2*M_PI);
                // Compute the gradient magnitude and apply a Gaussian weighing to give less emphasis to distant sample
                double t = lambda_descr*sigma_key;
                double M;
                if(gauss_flag)
                    M = hypot(dX, dY) * exp(-(X*X+Y*Y)/(2*t*t));
                else
                    M = hypot(dX, dY);

                // Determine the histogram and bin indices, Compute the (tri)linear weightings ...
                _myfloat alpha = X/(2*lambda_descr*sigma_key/Nhist) + (Nhist-1.0)/2.0;
                _myfloat beta  = Y/(2*lambda_descr*sigma_key/Nhist) + (Nhist-1.0)/2.0;
                _myfloat gamma = ori/(2*M_PI)*Nbins;
                //    ...and add contributions to respective bins in different histograms.
                // a loop with 1 or two elements
                int i0 = floor(alpha);
                int j0 = floor(beta);
                for(int i = MAX(0,i0);i<=MIN(i0+1,Nhist-1);i++){
                    for(int j = MAX(0,j0);j<=MIN(j0+1,Nhist-1);j++){ // looping through all surrounding histograms.

                        int k;
                        // Contribution to left bin.
                        k = ((int)gamma+Nbins)%Nbins;
                        descr[i*Nhist*Nbins+j*Nbins+k] += (1.-(gamma-floor(gamma)))
                                                         *(1.0-ABS((float)i-alpha))
                                                         *(1.0-ABS((float)j-beta))
                                                         *M;

                        // Contribution to right bin.
                        k = ((int)gamma+1+Nbins)%Nbins;
                        descr[i*Nhist*Nbins+j*Nbins+k] += (1.0-(floor(gamma)+1-gamma))
                                                         *(1.0-ABS((float)i-alpha))
                                                         *(1.0-ABS((float)j-beta))
                                                         *M;
                    }
                }
            }
        }
    }
}


static void ellipse_keypoints_attribute_oriented_descriptors(struct sift_scalespace *sx,
                                                             struct sift_scalespace *sy, /*gradient scalespaces*/
                                                             struct sift_keypoints *keys,
                                                             int n_bins,
                                                             int n_hist,
                                                             int n_ori,
                                                             _myfloat lambda_ori,
                                                             _myfloat lambda_descr,
                                                             int gauss_flag)
{
    for(int k = 0; k < keys->size; k++){

        /** Loading keypoint gradient scalespaces */
        struct keypoint* key = keys->list[k];
        _myfloat x = key->x;
        _myfloat y = key->y;
        int o = key->o;
        int s = key->s;
        _myfloat a = key->a;
        _myfloat b = key->b;
        _myfloat c = key->c;

        // load scalespace gradient
        int w = sx->octaves[o]->w;
        int h = sx->octaves[o]->h;
        _myfloat delta  = sx->octaves[o]->delta;
        const _myfloat* dx = &(sx->octaves[o]->imStack[s*w*h]);
        const _myfloat* dy = &(sy->octaves[o]->imStack[s*w*h]);

        //conversion to the octave's coordinates
        x /= delta;
        y /= delta;
        a *= (delta*delta);
        b *= (delta*delta);
        c *= (delta*delta);

        /** Accumulate gradient orientation histogram */
        ellipse_sift_accumulate_orientation_histogram(x, y,
                                                      a, b, c,
                                                      dx, dy, w, h,
                                                      n_bins,lambda_ori,key->orihist);

        /** Extract principal orientation (includes histogram smoothing)*/
         key->theta = sift_extract_one_orientation(key->orihist, key->n_bins);

        _myfloat theta = key->theta;
        /** Compute descriptor representation */
        ellipse_sift_extract_feature_vector(x, y,
                                            a, b, c,
                                            theta,
                                            dx, dy, w, h,
                                            n_hist, n_ori, lambda_descr, gauss_flag, key->descr);  /* output feature vector */

        /** Threshold and quantization of the descriptor */
        int n_descr = n_hist*n_hist*n_ori;
        sift_threshold_and_quantize_feature_vector(key->descr, n_descr, 0.2);
    }
}


void ellipse_sift_anatomy_orientation_and_description(const _myfloat* x,
                                              int w,
                                              int h,
                                              const struct sift_parameters* p,
                                              struct sift_keypoints* k,
                                              int gauss_flag) // TODO maybe include in parameter p
{
    // number octaves, limited by the maximum number of subsamplings
    int n_oct = MIN(p->n_oct, (int)(log(MIN(w,h)/p->delta_min)/M_LN2));

    // MEMORY ALLOCATION
    struct sift_scalespace* s   = sift_malloc_scalespace_lowe(n_oct,p->n_spo,w,h,p->delta_min,p->sigma_min);
    struct sift_scalespace* sx  = sift_malloc_scalespace_from_model(s);
    struct sift_scalespace* sy  = sift_malloc_scalespace_from_model(s);
    scalespace_compute(s, x, w, h, p->sigma_in);  // Builds Lowe's scale-space
    scalespace_compute_gradient(s,sx,sy); // Pre-computes gradient scale-space

    ellipse_keypoints_attribute_oriented_descriptors(sx, sy, k, p->n_bins, p->n_hist,p->n_ori, p->lambda_ori, p->lambda_descr, gauss_flag);


    sift_free_scalespace(s);
    sift_free_scalespace(sx);
    sift_free_scalespace(sy);
}


/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/



// Read a list of ellipses
//  WARNING - Convention changes (MIKO -> LOWE)
void sift_read_ellipses(struct sift_keypoints* keys,
                         const char* name,
                         int n_hist,
                         int n_ori,
                         int n_bins,
                         int flag)
{
    size_t buffer_size = 1024 * 1024;  // 1MB buffer fir long lines.
    char* buffer = xmalloc(buffer_size); 
    FILE* stream = fopen(name,"r");
    char * r;
    r = fgets(buffer, buffer_size, stream); // skip two lines
    r = fgets(buffer, buffer_size, stream);
    (void)r;
    while(fgets(buffer, buffer_size, stream) != NULL){
        int pos = 0;
        int read = 0;
        struct keypoint* key = sift_malloc_keypoint(n_ori, n_hist, n_bins);
        // read coordinates
        sscanf(buffer+pos,"%f %f %f %f %f %n", &key->y // WARNING different convention
                                             , &key->x
                                             , &key->c // WARNING different convention
                                             , &key->b
                                             , &key->a // WARNING different convention
                                               //   , &key->theta
                                             , &read);
        pos+=read;
        if (flag > 0){
            // read descriptor
            for(int i = 0; i < n_hist*n_hist*n_ori; i++){
                sscanf(buffer+pos, "%f %n",&(key->descr[i]),&read);
                pos +=read;
            }
        }
        if (flag > 1){
            // read keypoint orientation histogram
            for(int i=0;i<n_bins;i++){
                sscanf(buffer+pos, "%f %n",&(key->orihist[i]),&read);
                pos += read;
            }
        }
        // add the keypoint to the list
        sift_add_keypoint_to_list(key,keys);
    }
    xfree(buffer);
}


void fprintf_one_ellipse(FILE* f, const struct keypoint* k, int n_descr, int n_bins, int flag)
{
    // coordinates
    // WARNING different convention
    fprintf(f,"%f %f %f %f %f", k->y, k->x, k->c, k->b, k->a);
    if (flag>0){
        // descriptor
        for(int n = 0; n < n_descr; n++){
            fprintf(f,"%i ", (int)k->descr[n]);
        }
    }
    if (flag>1){
        // orientation histogram
        for(int n = 0; n < n_bins; n++){
            fprintf(f,"%f ", k->orihist[n]);
        }
    }
}


void fprintf_ellipses(FILE* f, const struct sift_keypoints* keys, int flag)
{
    int n_hist = keys->list[0]->n_hist;
    int n_ori  = keys->list[0]->n_ori;
    int n_descr = n_hist*n_hist*n_ori;
    int n_bins = keys->list[0]->n_bins;
    int n = keys->size;
    fprintf(f,"%i \n%i \n", n_descr, n);
    for(int k = 0; k < n; k++){
        struct keypoint* key = keys->list[k];
        fprintf_one_ellipse(f, key, n_descr, n_bins, flag);
        fprintf(f,"\n");
    }
}


void sift_save_ellipses(const struct sift_keypoints* keys, const char* name, int flag)
{
    FILE* f = fopen(name,"w");
    if (!f){
        fatal_error("Failed to open %s for writing\n", name); // TODO check fprintf security
    }
    fprintf_ellipses(f, keys, flag);
    fclose(f);
}


void sift_print_ellipses(const struct sift_keypoints* keys, int flag)
{
    fprintf_ellipses(stdout, keys, flag);
}


void print_pairs_ellipses(const struct sift_keypoints *k1,
                          const struct sift_keypoints *k2)
{
    int n = k1->size;
    for(int i = 0; i < n ;i++){
        fprintf_one_ellipse(stdout, k1->list[i], 0, 0, 0);
        fprintf(stdout, " ");
        fprintf_one_ellipse(stdout, k2->list[i], 0, 0, 0);
        fprintf(stdout, "\n");
    }
}



void estimate_scalespace_index(struct sift_keypoints *l, struct sift_parameters* p)
{
    int nspo = p->n_spo;
    _myfloat sigma_min = p->sigma_min;
    for(int k = 0; k < l->size; k++){


        struct keypoint* key = l->list[k];
        _myfloat a = key->a;
        _myfloat b = key->b;
        _myfloat c = key->c;
        _myfloat r1, r2, alpha;
        ellipse_cartesian_to_polar(a, b, c, &alpha, &r1, &r2);
        _myfloat sigma = sqrt(r1*r2);
        fprintf(stderr, "a %f b %f c %f // sigma_min %f / sigma %f \n  ", a, b, c, sigma_min, sigma); 


        int o,s;
        int n = (int)(round( nspo * log(sigma / sigma_min) /M_LN2  ));
        o = (n-1)/nspo;
        if (o > -1){
            s = (n-1)%nspo + 1;
        }
        else{
            o = 0;
            s = 0;
        }
        printf( "octave o %i and scale s %i\n", o, s); 
        key->o = o;
        key->s = s;
        key->sigma = sigma; // pas necessaire
    }
}






