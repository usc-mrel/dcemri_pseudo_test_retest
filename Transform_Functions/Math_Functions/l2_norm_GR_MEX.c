#include "mex.h"
#include <math.h>

#ifdef __GNU__
    #include <omp.h>
#endif

#ifndef MAXCORES
  #define MAXCORES 1
#endif 

void mexFunction(int nlhs, mxArray *left[], int nrhs, const mxArray *right[]) {
  
  /*  Declare variables */
  mwSize  elem;
  long long i;
  mxClassID precision,precision1;
  const mwSize size[]={1,1};
  mxArray *X1, *X2, *T, *Y;
  double  *pX1r, *pX1i, *pX2r, *pX2i, *pYr, *pT, Td;
  double  xr, xi, L2 = 0.0;
	float   *pX1rf, *pX1if, *pX2rf, *pX2if, *pYrf, *pTf, Tf;
  float   xrf, xif;
  
  /*  Get number of elements */
  elem = mxGetNumberOfElements(right[0]);
  
	/*  Obtain class */
	precision = mxGetClassID(right[0]);
  precision1 = mxGetClassID(right[1]);
	
  /*  Throw error if complexities mismatch */
  if (mxIsComplex(right[0]) != mxIsComplex(right[1]))
    mexErrMsgTxt("l1_norm_GR_MEX: Inputs real/complex mismatch");
  if (precision != precision1)
    mexErrMsgTxt("l1_norm_GR_MEX: Input data type mismatch");
  
	/*  Create output matrix */
	Y = mxCreateNumericArray(2, size, precision, mxREAL);
	
  /*  Get pointers to input and output arrays */
	if (precision == mxDOUBLE_CLASS) {
		pX1r = mxGetPr(right[0]);
		pX2r = mxGetPr(right[1]);
		pX1i = mxGetPi(right[0]);
		pX2i = mxGetPi(right[1]);
		pYr  = mxGetPr(Y);
	}
	else {
		pX1rf = mxGetData(right[0]);
		pX2rf = mxGetData(right[1]);
		pX1if = mxGetImagData(right[0]);
		pX2if = mxGetImagData(right[1]);
		pYrf  = mxGetData(Y);
	}
  
  /*  Get pointer to input scalar */
  if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
    pT = mxGetData(right[2]);
  else
    pTf = mxGetData(right[2]);
  
  /*  Convert scale factor to correct data type */
  if (precision == mxDOUBLE_CLASS) {
    if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
      Td = pT[0];
    else
      Td = (double) pTf[0];
  }
  else {
    if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
      Tf = (float) pT[0];
    else
      Tf = pTf[0];
  }
  
#ifdef __GNU__
    /* Set number of threads */
    omp_set_num_threads(MAXCORES);
#endif
    
  /*  Loop through and compute the l2 norm of the combined coefficients */
	if (precision == mxDOUBLE_CLASS) {
		#pragma omp parallel for private(i,xr,xi) reduction(+: L2)
		for (i=0; i<elem; i++) {
			xr = pX1r[i] + Td*pX2r[i];
			xi = pX1i[i] + Td*pX2i[i];
			L2 += xr*xr + xi*xi;
		}
		pYr[0] = L2;
	}
	else {
		#pragma omp parallel for private(i,xrf,xif) reduction(+: L2)
		for (i=0; i<elem; i++) {
			xrf = pX1rf[i] + Tf*pX2rf[i];
			xif = pX1if[i] + Tf*pX2if[i];
			L2 += xrf*xrf + xif*xif;
		}
		pYrf[0] = L2;
	}
	
  /* Return values */
  left[0] = Y;
}
