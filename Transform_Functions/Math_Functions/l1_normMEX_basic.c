#include "mex.h"
#include <omp.h>
#include <math.h>
#include "emmintrin.h"
#include "xmmintrin.h"

void mexFunction(int nlhs, mxArray *left[], int nrhs, const mxArray *right[]) {
    
    /*  Declare variables */
    mwSize  elem, cmplx, cmplx1, cmplx2, cmplx3;
    long long i, elem2;
	const   mwSize size[]={1,1};
	mxClassID precision, precision1;
    mxArray *X1, *X2, *T, *Y;
    double  *pX1r, *pX1i, *pX2r, *pX2i, *pYr, *pT, *pSd, Td, Sd;
    double  xr, xi, L1;
	float   *pX1rf, *pX1if, *pX2rf, *pX2if, *pYrf, *pTf, *pSf, Tf, Sf;
	float   xrf, xif, L1f;
    
    /*  Get number of elements */
    elem = mxGetNumberOfElements(right[0]);
	/*  mexPrintf("elem: %i\n",elem);*/
    
    /*  Test for complex and obtain data class */
    cmplx  = mxIsComplex(right[0]);
    cmplx1 = mxIsComplex(right[1]);
    cmplx2 = mxIsComplex(right[2]);
    cmplx3 = mxIsComplex(right[3]);
    if (cmplx != cmplx1)
        mexErrMsgTxt("Inputs 0 and 1 have different complexity");
    if (cmplx2)
        mexErrMsgTxt("Input 2 is complex (must be real)");
    if (cmplx3)
        mexErrMsgTxt("Input 3 is complex (must be real)");
    
    /*  Obtain and test data class  */
	precision  = mxGetClassID(right[0]);
    precision1 = mxGetClassID(right[1]);
    if (precision != precision1)
        mexErrMsgTxt("Inputs 0 and 1 have different precision");
    
    /*  Get pointers to input arrays and create output array */
	Y = mxCreateNumericArray(2, size, precision, mxREAL);
	if (precision == mxDOUBLE_CLASS) {
		pX1r = mxGetPr(right[0]);
		pX2r = mxGetPr(right[1]);
		if (cmplx) {
			pX1i = mxGetPi(right[0]);
			pX2i = mxGetPi(right[1]);
		}
        pYr = mxGetPr(Y);
	}
	else {
		pX1rf = mxGetData(right[0]);
		pX2rf = mxGetData(right[1]);
		if (cmplx) {
			pX1if = mxGetImagData(right[0]);
			pX2if = mxGetImagData(right[1]);
		}
        pYrf = mxGetData(Y);
	}
	
    /*  Get pointer to input scalar */
    if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
        pT = mxGetData(right[2]);
    else
        pTf = mxGetData(right[2]);
    
    /*  Get pointer to smoothing factor */
	if (mxGetClassID(right[3]) == mxDOUBLE_CLASS)
        pSd  = mxGetData(right[3]);
    else
        pSf  = mxGetData(right[3]);
    
    /*  Convert scalars to same data type as input arrays */
    if (precision == mxDOUBLE_CLASS) {
        if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
            Td = (double)pT[0];
        else
            Td = (double)pTf[0];
        if (mxGetClassID(right[3]) == mxDOUBLE_CLASS)
            Sd = (double)pSd[0];
        else
            Sd = (double)pSf[0];
    }
    else {
        if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
            Tf = (float)pT[0];
        else
            Tf = (float)pTf[0];
        if (mxGetClassID(right[3]) == mxDOUBLE_CLASS)
            Sf = (float)pSd[0];
        else
            Sf = (float)pSf[0];
    }
    
    /*  Set number of threads */
    omp_set_num_threads(16);
    
    /*  Loop through and compute the abs of the combined coefficients then sum */
	if (precision == mxDOUBLE_CLASS) {
		if (cmplx) {
            #pragma omp parallel for private(i,xr,xi) reduction(+: L1)
			for (i=0; i<elem; i++) {
				xr = pX1r[i] + Td*pX2r[i];
				xi = pX1i[i] + Td*pX2i[i];
				L1 += sqrt(xr*xr + xi*xi + Sd);
			}
		}
		else {
            #pragma omp parallel for private(i,xr) reduction(+: L1)
			for (i=0; i<elem; i++) {
                xr = pX1r[i] + Td*pX2r[i];
                L1 += sqrt(xr*xr + Sd);
				/*L1 += fabs(pX1r[i] + Td*pX2r[i]);*/
			}
		}
        pYr[0] = L1;
	}
    
	else {
		if (cmplx) {
			#pragma omp parallel for private(i,xrf,xif) reduction(+: L1)
			for (i=0; i<elem; i++) {
				xrf = pX1rf[i] + Tf*pX2rf[i];
				xif = pX1if[i] + Tf*pX2if[i];
				L1 += sqrt(xrf*xrf + xif*xif + Sf);
			}
            pYrf[0] = L1;
		}
		else {
            #pragma omp parallel for private(i,xrf) reduction(+: L1)
			for (i=0; i<elem; i++) {
                xrf = pX1rf[i] + Tf*pX2rf[i];
                L1 += sqrt(xrf*xrf + Sf);
				/*L1 += fabs(pX1rf[i] + Tf*pX2rf[i]);*/
			}
            pYrf[0] = L1;
		}
	}
	
    
    /* Return values */
    left[0] = Y;
}
