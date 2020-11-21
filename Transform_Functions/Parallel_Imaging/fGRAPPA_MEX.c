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
    mwSize  nD, iNP, iNV, iNS, iNT, iNRS;
    long long iNR;
	mwSize  indK, indG, indKG;
    mwSize  indKt, indGt, indKGt;
	mwSize  s4, s3, s2, s1, s0;
	const mwSize *sz;
	mxClassID precision;
    mxArray *K, *G, *KG;
    double  *pKr, *pKi, *pGr, *pGi, *pKGr, *pKGi;
	float   *pKrf, *pKif, *pGrf, *pGif, *pKGrf, *pKGif;
    
    /*  Get sizes */
	nD = mxGetNumberOfDimensions(right[0]);
    sz = mxGetDimensions(right[0]);
	s4 = sz[4]*sz[2]*sz[1]*sz[0]; /* not an error! */
	s3 = sz[3]*sz[2]*sz[1]*sz[0];
	s2 = sz[2]*sz[1]*sz[0];
	s1 = sz[1]*sz[0];
	s0 = sz[0];
	
	/*  Obtain data class */
	precision = mxGetClassID(right[0]);
    
    /*  Throw error if complexities/precision mismatch */
    if (mxIsComplex(right[0]) != mxIsComplex(right[1]))
        mexErrMsgTxt("fGRAPPA_MEX: Inputs real/complex mismatch");
    if (mxGetClassID(right[0]) != mxGetClassID(right[1]))
        mexErrMsgTxt("fGRAPPA_MEX: Input data type mismatch");
    
    /*  Create output, get pointers to input/output */
	KG = mxCreateNumericArray(nD,sz,precision,mxCOMPLEX);
	if (precision == mxDOUBLE_CLASS) {
		pKr = mxGetPr(right[0]);
		pKi = mxGetPi(right[0]);
		pGr = mxGetPr(right[1]);
		pGi = mxGetPi(right[1]);
		pKGr = mxGetPr(KG);
		pKGi = mxGetPi(KG);
	}
	else {
		pKrf = mxGetData(right[0]);
		pKif = mxGetImagData(right[0]);
		pGrf = mxGetData(right[1]);
		pGif = mxGetImagData(right[1]);
		pKGrf = mxGetData(KG);
		pKGif = mxGetImagData(KG);
	}
    
    
    #ifdef __GNU__
        /* Set number of threads */
        omp_set_num_threads(MAXCORES);
	#endif
    
	
	/*	Loop through elements */
	if (precision == mxDOUBLE_CLASS) {
		for (iNRS=0; iNRS<sz[4]; iNRS++) {
			#pragma omp parallel for private(iNR,iNT,iNS,iNV,iNP,indK,indG,indKG,indKGt,indKt,indGt)
			for (iNR=0; iNR<sz[4]; iNR++) {
			for (iNT=0; iNT<sz[3]; iNT++) {
			for (iNS=0; iNS<sz[2]; iNS++) {
			for (iNV=0; iNV<sz[1]; iNV++) {
                indKGt = iNV*s0 + iNS*s1 + iNT*s2  + iNR*s3;
                indKt  = iNV*s0 + iNS*s1 + iNT*s2  + iNRS*s3;
                indGt  = iNV*s0 + iNS*s1 + iNRS*s2 + iNR*s4;
                for (iNP=0; iNP<sz[0]; iNP++) {
                    indKG = iNP + indKGt;
                    indK  = iNP + indKt;
                    indG  = iNP + indGt;
                    pKGr[indKG] += pKr[indK]*pGr[indG] - pKi[indK]*pGi[indG];
                    pKGi[indKG] += pKr[indK]*pGi[indG] + pKi[indK]*pGr[indG];
                }
			}
			}
			}
			}
		}
	}
	else {
		for (iNRS=0; iNRS<sz[4]; iNRS++) {
			#pragma omp parallel for private(iNR,iNT,iNS,iNV,iNP,indK,indG,indKG,indKGt,indKt,indGt)
			for (iNR=0; iNR<sz[4]; iNR++) {
			for (iNT=0; iNT<sz[3]; iNT++) {
			for (iNS=0; iNS<sz[2]; iNS++) {
			for (iNV=0; iNV<sz[1]; iNV++) {
                indKGt = iNV*s0 + iNS*s1 + iNT*s2  + iNR*s3;
                indKt  = iNV*s0 + iNS*s1 + iNT*s2  + iNRS*s3;
                indGt  = iNV*s0 + iNS*s1 + iNRS*s2 + iNR*s4;
                for (iNP=0; iNP<sz[0]; iNP++) {
                    indKG = iNP + indKGt;
                    indK  = iNP + indKt;
                    indG  = iNP + indGt;
                    pKGrf[indKG] += pKrf[indK]*pGrf[indG] - pKif[indK]*pGif[indG];
                    pKGif[indKG] += pKrf[indK]*pGif[indG] + pKif[indK]*pGrf[indG];
                }
			}
			}
			}
			}
		}
	}
	
	
	/* Take difference with input */
	if (precision == mxDOUBLE_CLASS) {
        #pragma omp parallel for private(iNR,iNT,iNS,iNV,iNP,indK)
        for (iNR=0; iNR<sz[4]; iNR++) {
            for (iNT=0; iNT<sz[3]; iNT++) {
            for (iNS=0; iNS<sz[2]; iNS++) {
            for (iNV=0; iNV<sz[1]; iNV++) {
            for (iNP=0; iNP<sz[0]; iNP++) {
                indK  = iNP + iNV*s0 + iNS*s1 + iNT*s2  + iNR*s3;
                pKGr[indK] -= pKr[indK];
                pKGi[indK] -= pKi[indK];
            }
            }
            }
            }
        }
	}
	else{
        #pragma omp parallel for private(iNR,iNT,iNS,iNV,iNP,indK)
        for (iNR=0; iNR<sz[4]; iNR++) {
            for (iNT=0; iNT<sz[3]; iNT++) {
            for (iNS=0; iNS<sz[2]; iNS++) {
            for (iNV=0; iNV<sz[1]; iNV++) {
            for (iNP=0; iNP<sz[0]; iNP++) {
                indK  = iNP + iNV*s0 + iNS*s1 + iNT*s2  + iNR*s3;
                pKGrf[indK] -= pKrf[indK];
                pKGif[indK] -= pKif[indK];
            }
            }
            }
            }
        }
	}
	
    /* Return values */
    left[0] = KG;
}
