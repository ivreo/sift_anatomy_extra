#include "bsplines.h"


double pow_int(double z, int exponent){
	double out=1;
	for(int i=1;i<=exponent;i++){
		out *= z;
	}
	return out;
}

double init_causal_filter(double *in, int length, double pole){
	double init;
	init = in[0] + in[length-1]*pow_int(pole, length-1); //unduplicated coefficients of the extension.
	for(int n=1; n <= length-2 ;n++){
		init += in[n]*( pow_int(pole,n) + pow_int(pole, 2*length-2-n) ); //duplicated coefficients.
	}
	init = init/(1-pow_int(pole, 2*length-2));
	return init;
}


double init_anticausal_filter(double *in, int length, double pole){
	return (-pole/1-pole*pole)*(pole*in[length-2]+in[length-1]);
}

//void bspline_decomposition_1D(double* in, int length, double* poles, int nbpoles){


void bspline_decomposition_1D(double *in, int length, double *poles, int nbpoles){
	
	int p;
	int i;
	double normFactor = 1;
	
	// Normalization factor
	for(p = 0; p < nbpoles; p++){
		normFactor *= (1.-poles[p])*(1.-1./poles[p]);
	}
	
	// Series of filters (loop on poles)
	
	for(p = 0; p < nbpoles; p++){
		/// forward recursion
		// initializing
		in[0]= init_causal_filter(in, length, poles[p]);
		// filtering
		for(i=1;i<length;i++){
			in[i] += poles[p]*in[i-1];
		}
		
		/// backward recursion
		// initializing
		in[length-1] = init_anticausal_filter(in,length,poles[p]);
		// filtering
		for(i=length-2;i>=0;i--){  //TEST moisan commance avec [i = length-1]
			//printf("DEBUUUUG BACKWARD index value i = %i, length %i \n",i,length);
			in[i] = poles[p]*( in[i+1] - in[i] ); 
		}
	}
	
	//Normalization
	for(i=0;i<length;i++){
		in[i] = normFactor*in[i];
	}	
		
		
}
	
///  ///////////////////////// IMAGE BSPLINE DECOMPOSITION ///////////////////////////////////////////

double* image_bspline_decomposition(double *in, int order, int width, int height){
	printf("INSIDE B SPLINE DECOMPOSITION 1\n");
	int i,j;
	double *imcopy = (double*)malloc(width*height*sizeof(double));      //      new double[width*height];  // local copy
	double *transp = (double*)malloc(width*height*sizeof(double));      //      new double[height*width]; // for transposition
	int nbpoles;
	double poles[5]; // that covers bsplines up to order 11.
  
	double* out = (double*)malloc(width*height*sizeof(double));      //  new double[width*height];
  
	// initialize poles
	switch(order){
		case(2) : poles[0] = -0.17157288; /* sqrt(8)-3 */ nbpoles = 1;
			break;
			
		case 3: poles[0]=-0.26794919;  /* sqrt(3)-2 */  nbpoles = 1;
			break;

		case 4: poles[0]=-0.361341; poles[1]=-0.0137254; nbpoles = 2;
			break;

		case 5: poles[0]=-0.430575; poles[1]=-0.0430963; nbpoles = 2;
			break;
      
		case 6: poles[0]=-0.488295; poles[1]=-0.0816793; poles[2]=-0.00141415; nbpoles = 3;
			break;

		case 7: poles[0]=-0.53528; poles[1]=-0.122555; poles[2]=-0.00914869; nbpoles = 3;
			break;
      
		case 8: poles[0]=-0.574687; poles[1]=-0.163035; poles[2]=-0.0236323; poles[3]=-0.000153821; nbpoles = 4;
			break;

		case 9: poles[0]=-0.607997; poles[1]=-0.201751; poles[2]=-0.0432226; poles[3]=-0.00212131; nbpoles = 4;
			break;
      
		case 10: poles[0]=-0.636551; poles[1]=-0.238183; poles[2]=-0.065727; poles[3]=-0.00752819; poles[4]=-0.0000169828; nbpoles = 5;
			break;
      
		case 11: poles[0]=-0.661266; poles[1]=-0.27218; poles[2]=-0.0897596; poles[3]=-0.0166696;  poles[4]=-0.000510558; nbpoles = 5;
			break;
      
		default:
			printf("finvspline: order should be in 2..11.\n");
			exit(-1);
	}
	
// 	printf("bsplines ICI 1\n");
	/// Copy input image.
	for(i=0;i<width*height;i++){
		imcopy[i] = in[i];
	}
	
	
// 	printf("bsplines ICI 2\n");
	/// Decomposition on each line.
	for (i=0;i<height;i++){
		bspline_decomposition_1D(&imcopy[i*width], width, poles, nbpoles);
	}
	
	/// Transpose
	for (i=0;i<height;i++){
		for (j=0;j<width;j++){
			transp[j*height+i] = imcopy[i*width+j];			
		}
	}
// 	printf("bsplines ICI 3\n");
	/// Decomposition on each line (column)
	for (i=0;i<width;i++){
		bspline_decomposition_1D(&transp[i*height], height, poles, nbpoles);
	}	
// 	printf("bsplines ICI 4\n");
	/// Transpose
	for (i=0;i<height;i++){
		for (j=0;j<width;j++){
			out[i*width+j] = transp[j*height+i];
		}
	}
// 	printf("bsplines ICI 5\n");
	/// Freeing memory 
	free(transp);
	free(imcopy);
	
	return out;
}







/// ///////////////////// IMAGE RECONSTRUCTION ////////////////////////////////

double* precompute_bspline_binomials(int order){
	
	double* binomials =   (double*)malloc((2*(int)((order+1)/2)+2)*sizeof(double));      //     double[ 2*(int)((order+1)/2) +2 ];
	binomials[0] = 1;
	int k;
	for(k=0; k <= 2*(int)((order+1)/2); k++){
		binomials[k+1] = - (order-k)/(k+1)*binomials[k];
	}
	
	double tmp_order = order;
	double factorial=1;
	while(tmp_order>1){
		factorial *=tmp_order;
		tmp_order--;
	}
	for(k=0; k <= 2*(int)((order+1)/2)+1; k++)
		binomials[k] /=factorial;
	
	return binomials;
}

/// WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
/// BUG (using the explicit formula instead for the moment.)
double bspline_reconstruction_1D(double x, int order, double* binomials, double* bspline_decomp_1D){
	
	double bn;
	double ux;
	int k,kp;
	
	ux=0;
	for( kp = - (order+1)/2 ; kp <=  (order+1)/2 +1 ; kp++){
		bn=0;
		for (k=kp; k<= (order+1)/2 +1 ; k++){
			bn += binomials[k-kp] * pow( x-(int)x - k + (double)(order+1)/2.0 , order);
		}
		ux += bn*bspline_decomp_1D[kp+(int)x];
	}
	return ux;
}


/*

double  image_bspline_reconstruction(double x, double y, int order, double* binomials, double* bspline_decomp_image, int width, int height){
	
	(void)height;
	double bnx, bny;
	double res;
	int k,kp;
	int l,lp;
	
	res=0;
	for( kp = - (order+1)/2 ; kp <=  (order+1)/2 +1 ; kp++){
		bnx = 0;
		for (k=kp; k<= (order+1)/2 +1 ; k++){
			bnx += binomials[k-kp] * pow( x-(int)x - k + (double)(order+1)/2.0 , order);
		}
		for( lp = - (order+1)/2 ; lp <=  (order+1)/2 +1 ; lp++){
			bny = 0;
			for (l=lp; l<= (order+1)/2 +1 ; l++){
				bny += binomials[l-lp] * pow( y-(int)y - l + (double)(order+1)/2.0 , order);
			}
			res += bspline_decomp_image[kp*width+lp]*bnx*bny;
		}
	}
	return res;
}
*/






/// SOLUTION TEMPORAIRE RECONSTRUCTION B-SPLINES code GETREUER , ON FIXE UN ORDRE  7 ET ON REFORMULE L'IMPLEMENTATION DE LA RECONSTRUCTION UN AUTRE JOUR.








/** @brief Septic B-spline kernel (KernelRadius = 4) */
double BSpline7Kernel(double x){
    x = fabs(x);

    if(x <= 1)
    {
        double xSqr = x*x;
        return ((((35*x - 140)*xSqr + 560)*xSqr - 1680)*xSqr + 2416) / 5040;
    }
    else if(x <= 2)
    {
        x = 2 - x;
        return (120 + (392 + (504 + (280 + (-84 + (-42 +
            21*x)*x)*x*x)*x)*x)*x) / 5040;
    }
    else if(x < 3)
    {
        x = 3 - x;
        return (((((((-7*x + 7)*x + 21)*x + 35)*x + 35)*x
            + 21)*x + 7)*x + 1) / 5040;
    }
    else if(x < 4)
    {
        double xSqr;
        x = 4 - x;
        xSqr = x*x;
        return xSqr*xSqr*xSqr*x / 5040;
    }
    else
        return 0;
}

/** @brief Cubic B-spline kernel (KernelRadius = 2) */
double BSpline3Kernel(double x){
    x = fabs(x);

    if(x < 1)
        return (x/2 - 1)*x*x + 0.66666666666666667f;
    else if(x < 2)
    {
        x = 2 - x;
        return x*x*x/6;
    }
    else
        return 0;
}


// Enorme triffoullis 
double  image_bspline_reconstruction(double x, double y, int order, double* binomials, double* bspline_decomp_image, int width, int height){

	int kp,lp;
	int i,j;
	
	double res=0;
	for( kp = - order ; kp <=  order ; kp++){ // pas bon mais de toute façon les fonctions BSpline3Kernel() gerent le support
		i = kp+(int)x;
		while(i<0){i=i+(2*height-2);}
		while(i>2*height-2){i = i-(2*height-2);}
		if(i>height-1){i = (2*height-2)-i;}
		for( lp = - order ; lp <=  order ; lp++){
			//i = kp+(int)x; //printf("i %i \n",i);
			j = lp+(int)y; //printf("j %i \n",j);
			while(j<0){
                j=j+(2*width-2);
            }
			while(j>2*width-2){j = j - (2*width-2);}
			if(j>width-1){j = (2*width-2)-j;}
					
			//res += bspline_decomp_image[i*width+j] * BSpline7Kernel(x-(int)x-kp) * BSpline7Kernel(y-(int)y-lp);
			res += bspline_decomp_image[i*width+j] * BSpline3Kernel(x-(int)x-kp) * BSpline3Kernel(y-(int)y-lp);
		}
	}
	return res;	
}




/// origine

// // // /*
// // // 
// // // 
// // // 
// // // /** @brief Quadratic B-spline kernel (KernelRadius = 1.5) */
// // // static float BSpline2Kernel(float x)
// // // {
// // //     x = fabs(x);
// // // 
// // //     if(x <= 0.5f)
// // //         return 0.75f - x*x;
// // //     else if(x < 1.5f)
// // //     {
// // //         x = 1.5f - x;
// // //         return x*x/2;
// // //     }
// // //     else
// // //         return 0;
// // // }
// // // 
// // // 
// // // 
// // // 
// // // /** @brief Cubic B-spline kernel (KernelRadius = 2) */
// // // static float BSpline3Kernel(float x)
// // // {
// // //     x = fabs(x);
// // // 
// // //     if(x < 1)
// // //         return (x/2 - 1)*x*x + 0.66666666666666667f;
// // //     else if(x < 2)
// // //     {
// // //         x = 2 - x;
// // //         return x*x*x/6;
// // //     }
// // //     else
// // //         return 0;
// // // }
// // // 
// // // 
// // // 
// // // /** @brief Quintic B-spline kernel (KernelRadius = 3) */
// // // static float BSpline5Kernel(float x)
// // // {
// // //     x = fabs(x);
// // // 
// // //     if(x <= 1)
// // //     {
// // //         float xSqr = x*x;
// // //         return (((-10*x + 30)*xSqr - 60)*xSqr + 66) / 120;
// // //     }
// // //     else if(x < 2)
// // //     {
// // //         x = 2 - x;
// // //         return (1 + (5 + (10 + (10 + (5 - 5*x)*x)*x)*x)*x) / 120;
// // //     }
// // //     else if(x < 3)
// // //     {
// // //         float xSqr;
// // //         x = 3 - x;
// // //         xSqr = x*x;
// // //         return xSqr*xSqr*x / 120;
// // //     }
// // //     else
// // //         return 0;
// // // }
// // // 
// // // 
// // // 
// // // 
// // // 
// // // 
// // // 
// // // 
// // // 
// // // 
// // // */
