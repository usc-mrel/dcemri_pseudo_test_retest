#include "mex.h"
#include <math.h>
#include "emmintrin.h"
#include "xmmintrin.h"

#ifdef __GNU__
    #include <omp.h>
#endif


#ifndef MAXCORES
  #define MAXCORES 1
#endif 

void mexFunction(int nlhs, mxArray *left[], int nrhs, const mxArray *right[]) {
    
    /*  Declare variables */
    mwSize  elem, cmplx, cmplx1, cmplx2, cmplx3;
    long long i, elem2;
	const   mwSize size[]={1,1};
	mxClassID precision, precision1;
    mxArray *X1, *X2, *T, *Y;
    double  *pX1r, *pX1i, *pX2r, *pX2i, *pYr, *pT, *pSd, Td, Sd;
    double  xr, xi, L1, dL1[2];
    __m128d vTd, vSd, vL1, vxr, vxi;
	float   *pX1rf, *pX1if, *pX2rf, *pX2if, *pYrf, *pTf, *pSf, Tf, Sf;
	float   xrf, xif, L1f, dL1f[4];
    __m128  vTf, vSf, vL1f, vxrf, vxif;
    
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
    
    #ifdef __GNU__
        /* Set number of threads */
        omp_set_num_threads(MAXCORES);
	#endif
    
    /*  Loop through and compute the abs of the combined coefficients then sum */
	if (precision == mxDOUBLE_CLASS) {
		if (cmplx) {
			
            /*  Compute the number of elements for SIMD loop */
            elem2 = (elem/2)*2;
            
            /*  SIMD variables */
            vTd = _mm_load1_pd(&Td);
            vSd = _mm_load1_pd(&Sd);
            vL1 = _mm_setzero_pd();
            
            #pragma omp parallel for private(i,vxr,vxi) reduction(+: vL1)
            for (i=0; i<elem2; i+=2) {
                vxr = _mm_add_pd(_mm_load_pd(pX1r+i),_mm_mul_pd(vTd,_mm_load_pd(pX2r+i)));
                vxr = _mm_mul_pd(vxr,vxr);
                vxi = _mm_add_pd(_mm_load_pd(pX1i+i),_mm_mul_pd(vTd,_mm_load_pd(pX2i+i)));
                vxi = _mm_mul_pd(vxi,vxi);
                vL1 = _mm_add_pd(vL1,_mm_sqrt_pd(_mm_add_pd(_mm_add_pd(vxr,vxi),vSd)));
            }
            
            /*  Save results */
            _mm_store_pd(dL1,vL1);
            L1 = dL1[0] + dL1[1];
            
            /*  Finish the last few elements */
			for (i=elem2; i<elem; i++) {
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
            
            /*  Compute the number of elements for SIMD loop */
            elem2 = (elem/4)*4;
            
            /*  SIMD variables */
            vTf  = _mm_load1_ps(&Tf);
            vSf  = _mm_load1_ps(&Sf);
            vL1f = _mm_setzero_ps();
            
			#pragma omp parallel for private(i,vxrf,vxif) reduction(+: vL1f)
            for (i=0; i<elem2; i+=4) {
				vxrf = _mm_add_ps(_mm_load_ps(pX1rf+i),_mm_mul_ps(vTf,_mm_load_ps(pX2rf+i)));
                vxrf = _mm_mul_ps(vxrf,vxrf);
				vxif = _mm_add_ps(_mm_load_ps(pX1if+i),_mm_mul_ps(vTf,_mm_load_ps(pX2if+i)));
                vxif = _mm_mul_ps(vxif,vxif);
				vL1f = _mm_add_ps(vL1f,_mm_sqrt_ps(_mm_add_ps(_mm_add_ps(vxrf,vxif),vSf)));
			}
            
            /*  Save results */
            _mm_store_ps(dL1f,vL1f);
            L1f = dL1f[0] + dL1f[1] + dL1f[2] + dL1f[3];
            
            /*  Finish the last few elements */
			for (i=elem2; i<elem; i++) {
				xrf = pX1rf[i] + Tf*pX2rf[i];
				xif = pX1if[i] + Tf*pX2if[i];
				L1f += sqrt(xrf*xrf + xif*xif + Sf);
			}
            pYrf[0] = L1f;
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
