#include "mex.h"
#include <math.h>
#include "matrix.h"
#include <string.h>
#include "fast_mxArray_setup.c"

#ifdef __GNU__
    #include <omp.h>
#endif

#ifndef MAXCORES
  #define MAXCORES 1
#endif 

void mexFunction(int nlhs, mxArray *left[], int nrhs, const mxArray *right[]) {
    
    /*  Declare variables */
    mwSize nD, elem, cmplxX, cmplxY, *size2;
    long long i;
	const mwSize *size;
    mxComplexity cmplx;
	mxClassID precision1, precision3;
    mxArray *a, *X, *b, *Y, *Z;
    double  *pa, *pXr, *pXi, *pb, *pYr, *pYi, *pZr, *pZi, Ad=0.0, Bd=0.0;
	float   *paf, *pXrf, *pXif, *pbf, *pYrf, *pYif, *pZrf, *pZif, Af=0.0, Bf=0.0;
    
    /*  Get number of elements */
    elem = mxGetNumberOfElements(right[1]);
	nD   = mxGetNumberOfDimensions(right[1]);
	size = mxGetDimensions(right[1]);
    
    /*  Perform strange memory copy to replicate the size (needed for create_array_d/f) */
    size2 = (mwSize *)mxMalloc(nD*sizeof(mwSize));
    memcpy(size2,size,nD*sizeof(mwSize));
    
    /*  Test for complex arrays and obtain data class of input arrays */
    cmplxX = mxIsComplex(right[1]);
	cmplxY = mxIsComplex(right[3]);
	precision1 = mxGetClassID(right[1]);
    precision3 = mxGetClassID(right[3]);
    cmplx = cmplxX ? mxCOMPLEX:mxREAL;
    
    /*  Throw error if precision or complex mismatch */
    if (precision1 != precision3)
        mexErrMsgTxt("mult_add_MEX: Input data type mismatch");
    if (cmplxX != cmplxY)
        mexErrMsgTxt("Inputs 1 and 3 have different complexity");
    if (mxIsComplex(right[0]) | mxIsComplex(right[2]))
        mexErrMsgTxt("Input 0 and/or 2 is complex");
    
    
	/*  Create output array and pointers */
    if (precision1 == mxDOUBLE_CLASS)
        create_array_d(&(left[0]), &pZr, &pZi, nD, size2, cmplx, 0);
    else
        create_array_f(&(left[0]), &pZrf, &pZif, nD, size2, cmplx, 0);
    
    
    /*  Get pointers to input arrays */
	if (precision1 == mxDOUBLE_CLASS) {
		pXr = mxGetData(right[1]);
		pYr = mxGetData(right[3]);
		if (cmplxX) {
			pXi = mxGetImagData(right[1]);
			pYi = mxGetImagData(right[3]);
        }
	}
	else {
		pXrf = mxGetData(right[1]);
		pYrf = mxGetData(right[3]);
		if (cmplxX) {
			pXif = mxGetImagData(right[1]);
			pYif = mxGetImagData(right[3]);
        }
	}
    
    /*  Get pointers to input constants */
    if (mxGetClassID(right[0]) == mxDOUBLE_CLASS)
        pa  = mxGetData(right[0]);
    else
        paf = mxGetData(right[0]);
    if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
		pb  = mxGetData(right[2]);
    else
        pbf = mxGetData(right[2]);
    
    
	/*  Convert first scale factor to correct data type */
    if (precision1 == mxDOUBLE_CLASS) {
        /*mexPrintf("Converting A to double\n");*/
        if (mxGetClassID(right[0]) == mxDOUBLE_CLASS)
            Ad = pa[0];
        else
            Ad = (double) paf[0];
    }
    else {
        /*mexPrintf("Converting A to float\n");*/
        if (mxGetClassID(right[0]) == mxDOUBLE_CLASS)
            Af = (float) pa[0];
        else
            Af = paf[0];
    }
    /*  Convert second scale factor to correct data type */
    if (precision1 == mxDOUBLE_CLASS) {
        /*mexPrintf("Converting B to double\n");*/
        if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
            Bd = pb[0];
        else
            Bd = (double) pbf[0];
    }
    else {
        /*mexPrintf("Converting B to float\n");*/
        if (mxGetClassID(right[2]) == mxDOUBLE_CLASS)
            Bf = (float) pb[0];
        else
            Bf = pbf[0];
    }
    /*mexPrintf("Ad: %f\t Bd: %f\n",Ad,Bd);
     mexPrintf("Af: %f\t Bf: %f\n",(double) Af,(double) Bf);*/

    #ifdef __GNU__
        /* Set number of threads */
        omp_set_num_threads(MAXCORES);
	#endif
	
    /*  Loop through and compute the product */
	if (precision1 == mxDOUBLE_CLASS) {
		if (Ad == 1) {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZr[i] = pXr[i] + Bd*pYr[i];
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = pXi[i] + Bd*pYi[i];
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = pXi[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Bd*pYi[i];
			}
		}
		
		else if (Bd == 1) {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZr[i] = Ad*pXr[i] + pYr[i];
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Ad*pXi[i] + pYi[i];
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Ad*pXi[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = pYi[i];
			}
		}
		
		else if (Ad == -1) {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZr[i] = Bd*pYr[i] - pXr[i] ;
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Bd*pYi[i] - pXi[i] ;
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = -pXi[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Bd*pYi[i];
			}
		}
		
		else if (Bd == -1) {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZr[i] = Ad*pXr[i] - pYr[i];
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Ad*pXi[i] - pYi[i];
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Ad*pXi[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = -pYi[i];
			}
		}
		
		else {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZr[i] = Ad*pXr[i] + Bd*pYr[i];
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Ad*pXi[i] + Bd*pYi[i];
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Ad*pXi[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZi[i] = Bd*pYi[i];
			}
		}
	}
	
	
	else {
		if (Af == 1) {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZrf[i] = pXrf[i] + Bf*pYrf[i];
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = pXif[i] + Bf*pYif[i];
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = pXif[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Bf*pYif[i];
			}
		}
		
		else if (Bf == 1) {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZrf[i] = Af*pXrf[i] + pYrf[i];
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Af*pXif[i] + pYif[i];
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Af*pXif[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = pYif[i];
			}
		}
		
		else if (Af == -1) {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZrf[i] = Bf*pYrf[i] - pXrf[i] ;
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Bf*pYif[i] - pXif[i] ;
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = -pXif[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Bf*pYif[i];
			}
		}
		
		else if (Bf == -1) {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZrf[i] = Af*pXrf[i] - pYrf[i];
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Af*pXif[i] - pYif[i];
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Af*pXif[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = -pYif[i];
			}
		}
		
		else {
            #pragma omp parallel for private(i)
			for (i=0; i<elem; i++)
				pZrf[i] = Af*pXrf[i] + Bf*pYrf[i];
			if (cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Af*pXif[i] + Bf*pYif[i];
			}
			else if (cmplxX && !cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Af*pXif[i];
			}
			else if (!cmplxX && cmplxY) {
                #pragma omp parallel for private(i)
				for (i=0; i<elem; i++)
					pZif[i] = Bf*pYif[i];
			}
		}
	}
    
    /*  Free memory */
    mxFree(size2);
}
