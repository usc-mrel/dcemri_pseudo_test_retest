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
    mwSize  elem, cmpx0, cmpx1, cmpx2, cmpx3;
    long long i;
	const mwSize size[]={1,1};
	mxClassID precision0, precision1, precision2, precision3;
    mxArray *KN, *KD, *K, *T, *Y;
    double  *pKNr, *pKNi, *pKDr, *pKDi, *pKr, *pKi, *pY, *pT, Td;
    double  xr, xi, L2 = 0;
	float   *pKNrf, *pKNif, *pKDrf, *pKDif, *pKrf, *pKif, *pYf, *pTf, Tf;
    float   xrf, xif;
    
    /*  Get number of elements and obtain data class */
    elem = mxGetNumberOfElements(right[0]);
	precision0 = mxGetClassID(right[0]);
    precision1 = mxGetClassID(right[1]);
    precision2 = mxGetClassID(right[2]);
    
    /*  Throw error if precision mismatch */
    if ((precision0 != precision1) | (precision0 != precision2))
        mexErrMsgTxt("l1_norm_k_MEX: Input data type mismatch");
    
    /*  Test to ensure inputs are proper data type (real or complex) */
    if (!mxIsComplex(right[0]))
        mexErrMsgTxt("Input 0 (Knew) is real");
    if (!mxIsComplex(right[1]))
        mexErrMsgTxt("Input 1 (Kdelta) is real");
    if (!mxIsComplex(right[2]))
        mexErrMsgTxt("Input 2 (Kmeasured) is real");
    if (mxIsComplex(right[3]))
        mexErrMsgTxt("Input 3 (t) is complex");
	
    /*  Get pointers to input/output arrays and create output array */
	Y = mxCreateNumericArray(2, size, precision0, mxREAL);
	if (precision0 == mxDOUBLE_CLASS) {
		pKNr = mxGetPr(right[0]);
		pKNi = mxGetPi(right[0]);
		pKDr = mxGetPr(right[1]);
		pKDi = mxGetPi(right[1]);
		pKr  = mxGetPr(right[2]);
		pKi  = mxGetPi(right[2]);	
		pY   = mxGetPr(Y);
	}
	else {
		pKNrf = mxGetData(right[0]);
		pKNif = mxGetImagData(right[0]);
		pKDrf = mxGetData(right[1]);
		pKDif = mxGetImagData(right[1]);
		pKrf  = mxGetData(right[2]);
		pKif  = mxGetImagData(right[2]);	
		pYf   = mxGetData(Y);
	}
    
    /*  Get pointer to input scalar */
    if (mxGetClassID(right[3]) == mxDOUBLE_CLASS)
		pT  = mxGetData(right[3]);
    else
        pTf = mxGetData(right[3]);
    
    /*  Convert scale factor to correct data type */
    if (precision0 == mxDOUBLE_CLASS) {
        if (mxGetClassID(right[3]) == mxDOUBLE_CLASS)
            Td = pT[0];
        else
            Td = (double) pTf[0];
    }
    else {
        if (mxGetClassID(right[3]) == mxDOUBLE_CLASS)
            Tf = (float) pT[0];
        else
            Tf = pTf[0];
    }
    
    #ifdef __GNU__
        /* Set number of threads */
        omp_set_num_threads(MAXCORES);
	#endif
    
    /*  Loop through and compute the l2 norm of the combined coefficients */
	if (precision0 == mxDOUBLE_CLASS) {
        #pragma omp parallel for private(i,xr,xi) reduction(+: L2)
		for (i=0; i<elem; i++) {
			xr = (pKNr[i] + Td*pKDr[i] - pKr[i]);
			xi = (pKNi[i] + Td*pKDi[i] - pKi[i]);
			L2 += xr*xr + xi*xi;
		}
        pY[0] = L2;
	}
	else {
        #pragma omp parallel for private(i,xrf,xif) reduction(+: L2)
		for (i=0; i<elem; i++) {
			xrf = (pKNrf[i] + Tf*pKDrf[i] - pKrf[i]);
			xif = (pKNif[i] + Tf*pKDif[i] - pKif[i]);
			L2 += xrf*xrf + xif*xif;
		}
        pYf[0] = L2;
	}
    
    
    /* Return values */
    left[0] = Y;
}
