/**************************************************************************
 MEX function to compute the approximate gradient of the absolute value
 
 Author: R. Marc Lebel
 Contact: mlebel@gmail.com
 Date: 11/2010
 
 Useage: wc2 = absgradMEX(wc,smooth)
 
 Input:
 wc: numeric array (single/double; real/complex)
 smooth: small smoothing factor to prevent Inf
 
 Output:
 wc2: numeric array
**************************************************************************/

#include <stdio.h>
#include "mex.h"
#include <omp.h>
#include <math.h>
#include <string.h>
#include "fast_mxArray_setup.c"

float Q_rsqrt( float number )
{
  long i;
  float x2, y;
  const float threehalfs = 1.5F;
  
  x2 = number * 0.5F;
  y  = number;
  i  = * ( long * ) &y;
  i  = 0x5f375a86 - ( i >> 1 );
  y  = * ( float * ) &i;
  y  = y * ( threehalfs - ( x2 * y * y ) );
  y  = y * ( threehalfs - ( x2 * y * y ) );
  
  return y;
}

float Q_dsqrt( float number )
{
  long i;
  double x2, y;
  const double threehalfs = 1.5D;
  
  x2 = number * 0.5D;
  y  = number;
  i  = * ( long * ) &y;
  i  = 0x5fe6eb50c7b537a9 - ( i >> 1 );
  y  = * ( double * ) &i;
  y  = y * ( threehalfs - ( x2 * y * y ) );
  y  = y * ( threehalfs - ( x2 * y * y ) );
  y  = y * ( threehalfs - ( x2 * y * y ) );
  
  return y;
}


void mexFunction(int nlhs, mxArray *left[], int nrhs, const mxArray *right[]) {
    
    /*  Declare variables */
    mwSize  nD, elem, cmplx, *size2;
    long long i;
	mxClassID precision;
    const mwSize *size;
    mxComplexity cmplx2;
    mxArray *X, *Y;
    double  *pXr, *pXi, *pYi, *pYr, *pS, Sd, denom;
	float   *pXrf, *pXif, *pYif, *pYrf, *pSf, Sf, denomf;
    
    /*  Get size */
    nD   = mxGetNumberOfDimensions(right[0]);
    size = mxGetDimensions(right[0]);
    elem = mxGetNumberOfElements(right[0]);
    /*mexPrintf("nD: %i\n",nD);
     mexPrintf("size: %i\n",size[0]);
     mexPrintf("elem: %i\n",elem);*/
    
    /*  Perform strange memory copy to replicate the size (needed for create_array_d/f) */
    size2 = (mwSize *)mxMalloc(nD*sizeof(mwSize));
    memcpy(size2,size,nD*sizeof(mwSize));
    
    /*  Test for complex and obtain data class */
    cmplx = mxIsComplex(right[0]);
	precision = mxGetClassID(right[0]);
    cmplx2 = cmplx ? mxCOMPLEX:mxREAL;
    
    /*  Test to ensure smoothing factor is real */
    if (mxIsComplex(right[1]))
        mexErrMsgTxt("Inputs 1 is complex");
    
    /*  Get pointers to input array and create output */
	if (precision == mxDOUBLE_CLASS) {
		pXr = mxGetPr(right[0]);
		if (cmplx)
			pXi = mxGetPi(right[0]);
        
		/*  Create output and assign pointers */
        create_array_d(&(left[0]), &pYr, &pYi, nD, size2, cmplx2, 0);
	}
	else {
		pXrf = mxGetData(right[0]);
		if (cmplx)
			pXif = mxGetImagData(right[0]);
        
        /*  Create output and assign pointers */
        create_array_f(&(left[0]), &pYrf, &pYif, nD, size2, cmplx2, 0);
	}
    
    /*  Get pointer to input scalar */
	if (mxGetClassID(right[1]) == mxDOUBLE_CLASS)
        pS  = mxGetData(right[1]);
    else
        pSf  = mxGetData(right[1]);
    
    /*  Convert smoothing factor to appropriate class */
    if (precision == mxDOUBLE_CLASS) {
        if (mxGetClassID(right[1]) == mxDOUBLE_CLASS)
            Sd = pS[0];
        else
            Sd = (double) pSf[0];
    }
    else {
        if (mxGetClassID(right[1]) == mxDOUBLE_CLASS)
            Sf = (float) pS[0];
        else
            Sf = pSf[0];
    }
    
    
    
    /*  Loop through and compute the gradient of the absolute value */
	if (precision == mxDOUBLE_CLASS) {
		if (cmplx) {
            #pragma omp parallel for private(i,denom)
			for (i=0; i<elem; i++) {
                denom = 1.0/sqrt(pXr[i]*pXr[i] + pXi[i]*pXi[i] + Sd);
				pYr[i] = pXr[i] * denom;
				pYi[i] = pXi[i] * denom;
			}
		}
		else {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++) {
				pYr[i] = pXr[i]/sqrt(pXr[i]*pXr[i] + Sd);
			}
		}
	}
	else {
		if (cmplx) {
            #pragma omp parallel for private(i,denomf)
			for (i=0; i<elem; i++) {
                denomf = Q_rsqrt(pXrf[i]*pXrf[i] + pXif[i]*pXif[i] + Sf); /* Not working on ubuntu?! */
                /*denomf = 1.0/sqrt(pXrf[i]*pXrf[i] + pXif[i]*pXif[i] + Sf);*/
				pYrf[i] = pXrf[i] * denomf;
				pYif[i] = pXif[i] * denomf;
			}
		}
		else {
            #pragma omp parallel for private(i,denomf)
			for (i=0; i<elem; i++) {
				/*pYrf[i] = pXrf[i]/sqrt(pXrf[i]*pXrf[i] + Sf);*/
                denomf = Q_rsqrt(pXrf[i]*pXrf[i] + Sf);
                pYrf[i] = pXrf[i] * denomf; /* Not working on ubuntu?! */
			}
		}
	}
    
    
    /*  Free memory */
    mxFree(size2);
}
