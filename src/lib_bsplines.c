//#include "sp.h"
#include "lib_bsplines.h"
#include <math.h>   // for the declaration of pow
#include <string.h>  // for the declaration of memset (used here to initialize a vector to 0).



//%%%%%%%%%%%%%%%%%%%% GABRIELE's code

static _myfloat initcausal(_myfloat *c,int n,_myfloat z)
{
    _myfloat zk,z2k,iz,sum;
    int k;

    zk = z; iz = 1./z;
    z2k = pow(z,(_myfloat)n-1.);
    sum = c[0] + z2k * c[n-1];
    z2k = z2k*z2k*iz;
    for (k=1;k<=n-2;k++) {
        sum += (zk+z2k)*c[k];
        zk *= z;
        z2k *= iz;
    }
    return (sum/(1.-zk*zk));
}

_myfloat initanticausal(_myfloat *c,int n, _myfloat z)
{
    return((z/(z*z-1.))*(z*c[n-2]+c[n-1]));
}


static void invspline1D(_myfloat *c,int size,_myfloat *z,int npoles)
{
    _myfloat lambda;
    int n,k;

    /* normalization */
    for (k=npoles,lambda=1.;k--;) lambda *= (1.-z[k])*(1.-1./z[k]);
    for (n=size;n--;) c[n] *= lambda;

    /*----- Loop on poles -----*/
    for (k=0;k<npoles;k++) {

        /* forward recursion */
        c[0] = initcausal(c,size,z[k]);
        for (n=1;n<size;n++) 
            c[n] += z[k]*c[n-1];

        /* backwards recursion */
        c[size-1] = initanticausal(c,size,z[k]);
        for (n=size-1;n--;) 
            c[n] = z[k]*(c[n+1]-c[n]);

    }
}


/*------------------------------ MAIN MODULE ------------------------------*/
// Spline transform  (aka prefilter)
// in can be the same vector as out;
static void finvspline(const _myfloat *in, int nx, int order, _myfloat *out)
{
    _myfloat *c,*d,z[5];
    (void) d;
    (void) c;
    int npoles,x;

    //iv:  I removed pair polynomial order, It's not working properly

    /* initialize poles of associated z-filter */
    switch (order) 
    {
      //  case 2: z[0]=-0.17157288;  /* sqrt(8)-3 */
      //          break;

        //case 3: z[0]=-0.26794919;  /* sqrt(3)-2 */ 
        case 3: z[0]= -0.2679491924311228068233958765631541609764 ;  /* sqrt(3)-2 */ 
                break;

      //  case 4: z[0]=-0.361341; z[1]=-0.0137254;
      //          break;

        case 5: z[0]=-0.430575; z[1]=-0.0430963;
                break;

      //  case 6: z[0]=-0.488295; z[1]=-0.0816793; z[2]=-0.00141415;
      //          break;

        case 7: z[0]=-0.53528; z[1]=-0.122555; z[2]=-0.00914869;
                break;

      //  case 8: z[0]=-0.574687; z[1]=-0.163035; z[2]=-0.0236323; z[3]=-0.000153821;
      //          break;

        case 9: z[0]=-0.607997; z[1]=-0.201751; z[2]=-0.0432226; z[3]=-0.00212131;
                break;

      //  case 10: z[0]=-0.636551; z[1]=-0.238183; z[2]=-0.065727; z[3]=-0.00752819;
      //           z[4]=-0.0000169828;
      //           break;

        case 11: z[0]=-0.661266; z[1]=-0.27218; z[2]=-0.0897596; z[3]=-0.0166696; 
                 z[4]=-0.000510558;
                 break;

        default:
                 //fprintf(stderr,"finvspline: order should be in 2..11.\n");
                 fprintf(stderr,"finvspline: order should be in 3 5 7 9 11 \n");
    }
    npoles = order/2;

    /* initialize _myfloat array containing image */
    for (x=nx;x--;) 
        out[x] = (_myfloat)in[x];

    /* apply filter on lines */
    invspline1D(out,nx,z,npoles);

}



//    /*Linear interpolation*/
//    
//    static void linear(_myfloat *c,_myfloat t)
//    {
//        c[0]=t; c[1]=1.-t;
//    }
//    
//    /*Cubic Interpolation*/
//    
//    /* coefficients for cubic interpolant (Keys' function) */



//    static void keys(_myfloat *c,_myfloat t,_myfloat a)
//    {
//        _myfloat t2,at;
//    
//        t2 = t*t;
//        at = a*t;
//        c[0] = a*t2*(1.0-t);
//        c[1] = (2.0*a+3.0 - (a+2.0)*t)*t2 - at;
//        c[2] = ((a+2.0)*t - a-3.0)*t2 + 1.0;
//        c[3] = a*(t-2.0)*t2 + at;
//    }


/*BSpline-3 interpolation*/

/* coefficients for cubic spline */

static void spline3(_myfloat *c,_myfloat t)
{
    _myfloat tmp;

    tmp = 1.-t;
    c[0] = 0.1666666666*t*t*t;
    c[1] = 0.6666666666-0.5*tmp*tmp*(1.+t);
    c[2] = 0.6666666666-0.5*t*t*(2.-t);
    c[3] = 0.1666666666*tmp*tmp*tmp;

}


/** Indirect Interpolation  **/

/* pre-computation for spline of order >3 */
static void init_splinen(_myfloat *a,int n)

{

    a[0] = 1.;
    for (int k=2;k<=n;k++){ a[0]/=(_myfloat)k;}
    for (int k=1;k<=n+1;k++)
    { a[k] = - a[k-1] *(_myfloat)(n+2-k)/(_myfloat)k;}
}

_myfloat ipow(_myfloat x,int n)

{
    _myfloat res=0.0;
    res=pow(x,n);
    return res;
}

/* coefficients for spline of order >3 */
static void splinen(_myfloat *c,_myfloat t,_myfloat *a,int n)
{
    int i,k;
    _myfloat xn;

    memset((void *)c, 0, (n+1)*sizeof(_myfloat) );





    for (k=0;k<=n+1;k++) {
        xn = ipow(t+(_myfloat)k,n);
        for (i=k;i<=n;i++)
            c[i] += a[i-k]*xn;
    }
}


//         
//         // WTF !!
//         static int extansion(int N,int x)
//         {
//             if(x<0){x=-x;}
//             else if(x>=N){x=(2*N - 2) - x;}
//             return  x;
//         }
//         
//         





// //                      //  
// //                      //  
// //                      //  
/**************************************************************
 *
 *
 * MY STUFF (gabriele)
 *
 *
 * *********************/






/* extract image value (even outside image domain) */
_myfloat v(_myfloat *in, int n, int xi)
{
    // assume mirror boundaries 
    int x=xi;
    if(xi < 0)   x = -xi - 1;
    if(xi > n-1) x = -xi + 2*n - 1;

    return(in[x]);
}



_myfloat *SPLINE_EVAL_LUT =NULL;
int SPLINE_EVAL_LUT_STEPS = 0;

// TODO: receive a list of points
// TODO: FIXME errors at the boundaries..
void init_eval_LUT(int steps, int ord) {

    if (SPLINE_EVAL_LUT) free(SPLINE_EVAL_LUT);

    SPLINE_EVAL_LUT = (_myfloat*)malloc(12*steps*sizeof*SPLINE_EVAL_LUT);
    SPLINE_EVAL_LUT_STEPS = steps;

    for(int i=0;i<steps;i++) {
        _myfloat ak[13];
        _myfloat *c = SPLINE_EVAL_LUT + i*12;

        _myfloat u = (_myfloat)i/(_myfloat)steps;

        // evaluate spline to interpolate

        //int n1;
        int n2;
        if (ord==1){
            n2 = 1;
            c[0] = u;
            c[1] = 1. - u;
        }
        else if (ord==3) 
        {
            n2 = 2; 
            spline3(c,u); 
        } 
        else if (ord>3)
        {
            n2 = (1+ord)/2;
            init_splinen(ak,ord);
            splinen(c,u,ak,ord); 
        }
        // iv set but not used n2
        (void) n2;

    }

}



// TODO: receive a list of points
// TODO: FIXME errors at the boundaries..
_myfloat eval_spline(_myfloat * vec, int length, _myfloat xp, int ord) {
    _myfloat ak[13];
    _myfloat c[12];

    int  xi = (int)floor((_myfloat)xp);
    _myfloat u = xp-(_myfloat)xi;
    int n1,n2;

    if (ord==1){
        n2 = 1;
        c[0] = u;
        c[1] = 1. - u;
    }
    else if (ord==3) 
    {
        n2 = 2; 
        spline3(c,u); 
    } 
    else // ord>3
    {
        n2 = (1+ord)/2;
        init_splinen(ak,ord);  // TODO this is constant
        splinen(c,u,ak,ord); 
    }

    // Follows moisan's fcrop, where the cases near the boundaries are handled separately
    _myfloat res=0.;
    n1 = 1-n2;
    if (xi+n1>=0 && xi+n2<length) 
    {
        for (int d=n1;d<=n2;d++) {
            res += c[n2-d]*vec[xi+d];
        }
    }
    else
    {
        for (int d=n1;d<=n2;d++) 
            res += c[n2-d]*v(vec,length, xi+d);
    }
    return res;
}











// TODO: receive a list of points
// TODO: FIXME errors at the boundaries..
_myfloat eval_splineLUT(_myfloat * vec, int length, _myfloat xp, int ord) {

    int  xi = (int)floor((_myfloat)xp);
    _myfloat u = xp-(_myfloat)xi;
    int n1,n2;

    int idx = SPLINE_EVAL_LUT_STEPS*u;
    _myfloat *c = SPLINE_EVAL_LUT + (12*idx);

    if (ord==1) { n2 = 1; } 
    else if (ord==3) { n2 = 2; } 
    else { n2 = (1+ord)/2; }


    // Follows moisan's fcrop, where the cases near the boundaries are handled separately
    _myfloat res=0.;
    n1 = 1-n2;
    if (xi+n1>=0 && xi+n2<length) 
    {
        for (int d=n1;d<=n2;d++) 
            res += c[n2-d]*vec[xi+d];
    } 
    else  
    {
        for (int d=n1;d<=n2;d++) 
            res += c[n2-d]*v(vec,length, xi+d);
    }

    return res;
}




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END GABRIELE's code







void oversample_bsplines(const _myfloat* x, int wx, int hx,
                               _myfloat* y, int wy, int hy,
                               _myfloat delta, int order)
{
    // interpolation, horizontal direction
    _myfloat* tmp = xmalloc(hx*wy*sizeof(_myfloat));
    for( int i = 0; i < hx; i++){
        _myfloat* f = xmalloc(wx*sizeof(_myfloat));
        finvspline( &x[i*wx], wx , order, f);
        for( int j = 0; j < wy; j++){
            tmp[i*wy+j] = eval_spline(f, wx, delta*j, order);
        }
        free(f);
    }
    // interpolation, vertical direction
    for( int j = 0; j < wy; j++){
        _myfloat * tt = xmalloc(hx*sizeof(_myfloat));
        for(int i = 0; i < hx; i++){
            tt[i]  = tmp[i*wy+j];
        }
        _myfloat* f = xmalloc(wx*sizeof(_myfloat));
        finvspline( tt, hx , order, f);
        for( int i = 0; i < hy; i++){
            y[i*wy+j] = eval_spline(f, hx, delta*i, order);
        }
        free(f);
    }
}





//void translation_bsplines(_myfloat* x, _myfloat* y, int w, int h, double dx, double dy, int order)
//{
//    double* tmp = (double*)malloc(w*h*sizeof(double));
//    // dy
//    for(int i = 0; i < h; i++){
//        double* f = (double*)malloc(w*sizeof(double));
//        finvspline( &x[i*w], w , order, f);
//        for( int j = 0; j < w; j++){
//            tmp[i*w+j] = eval_spline(f, w, j+dy, order);
//        }
//        free(f);
//    }
//    // dx 
//    for( int j = 0; j < w; j++){
//        // transposition
//        double * tt = (double*)malloc(h*sizeof(double));
//        for(int i = 0; i < h; i++){
//            tt[i]  = tmp[i*w+j];
//        }
//        double* f = (double*)malloc(w*sizeof(double));
//        finvspline(tt, h , order, f);
//        for( int i = 0; i < h; i++){
//            y[i*w+j] = eval_spline(f, h, i + dx, order);
//        }
//        free(f);
//    }
//}
