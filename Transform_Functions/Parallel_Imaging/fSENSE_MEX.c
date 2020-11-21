/**************************************************************************
 MEX function to compute the approximate gradient of the absolute value
 
 Author: R. Marc Lebel
 Contact: mlebel@gmail.com
 Date: 11/2010
 
 Useage: imgS = absgradMEX(img,sens)
 
 Input:
 img: numeric array (single/double; real/complex)
 sens: sensitivity maps
 
 Output:
 imgS: multichanel image
**************************************************************************/

#include "mex.h"
#include <math.h>
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
    mwSize i, j, k, t;
    long long r;
    mwSize np, nv, ns, nt, nr;
    mwSize indS, indIMG, indIMGS, indXY, indXYZ, indXYZT, indT1, indT2, indT3;
    mwSize  nD, *sizeOUT;
	mxClassID precision;
    double  *pSZd, *pIMGrd, *pIMGid, *pSrd, *pSid, *pIMGSrd, *pIMGSid;
	float   *pSZf, *pIMGrf, *pIMGif, *pSrf, *pSif, *pIMGSrf, *pIMGSif;
    
    /*  Get size */
    nD = 5;
    /*mexPrintf("nD: %i\n",nD); */
    if (mxGetClassID(right[2]) == mxDOUBLE_CLASS) {
        pSZd = mxGetPr(right[2]);
        np = (int)pSZd[0]; nv = (int)pSZd[1]; ns = (int)pSZd[2]; nt = (int)pSZd[3]; nr = (int)pSZd[4];
    }
    else {
        pSZf = mxGetData(right[2]);
        np = (int)pSZf[0]; nv = (int)pSZf[1]; ns = (int)pSZf[2]; nt = (int)pSZf[3]; nr = (int)pSZf[4];
    }
    /*mexPrintf("size: %i x %i x %i x %i x %i\n",np,nv,ns,nt,nr);*/
    
    /*  Perform strange memory copy to replicate the size (needed for create_array_d/f) */
    sizeOUT = (mwSize *)mxMalloc(nD*sizeof(mwSize));
    sizeOUT[0] = np; sizeOUT[1] = nv; sizeOUT[2] = ns; sizeOUT[3] = nt; sizeOUT[4] = nr;
    /*mexPrintf("sizeOUT: %i x %i x %i x %i x %i\n",sizeOUT[0],sizeOUT[1],sizeOUT[2],sizeOUT[3],sizeOUT[4]);*/
    
    /*  Test for complex and obtain data class */
    if (!(mxIsComplex(right[0]) & mxIsComplex(right[1])))
        mexErrMsgTxt("Inputs need to be complex");
   	precision = mxGetClassID(right[0]);
    
    /*  Get pointers to input array and create output */
	if (precision == mxDOUBLE_CLASS) {
		pIMGrd = mxGetPr(right[0]);
        pIMGid = mxGetPi(right[0]);
        pSrd   = mxGetPr(right[1]);
        pSid   = mxGetPi(right[1]);
        
		/*  Create output and assign pointers */
        /*create_array_d(&(left[0]), &pIMGSrd, &pIMGSid, nD, sizeOUT, mxCOMPLEX, 0);*/
        left[0] = mxCreateNumericArray(nD,sizeOUT,precision,mxCOMPLEX);
        pIMGSrd = mxGetPr(left[0]);
        pIMGSid = mxGetPi(left[0]);
	}
	else {
		pIMGrf = mxGetData(right[0]);
        pIMGif = mxGetImagData(right[0]);
        pSrf   = mxGetData(right[1]);
        pSif   = mxGetImagData(right[1]);
        
        /*  Create output and assign pointers */
        /*create_array_f(&(left[0]), &pIMGSrf, &pIMGSif, nD, sizeOUT, mxCOMPLEX, 0);*/
        left[0] = mxCreateNumericArray(nD,sizeOUT,precision,mxCOMPLEX);
        pIMGSrf = mxGetData(left[0]);
        pIMGSif = mxGetImagData(left[0]);
	}
    
    #ifdef __GNU__
        /* Set number of threads */
        omp_set_num_threads(MAXCORES);
	#endif
    
    /*  Loop through elements */
    indXY   = np*nv;
    indXYZ  = np*nv*ns;
    indXYZT = np*nv*ns*nt;
    if (precision == mxDOUBLE_CLASS) {
        #pragma omp parallel for private(i,j,k,t,r,indT1,indT2,indT3,indIMG,indIMGS,indS)
        for (r=0; r<nr; r++) {
            for (t=0; t<nt; t++) {
                for (k=0; k<ns; k++) {
                for (j=0; j<nv; j++) {
                    indT1 =             t*indXYZ + k*indXY + j*np;
                    indT2 = r*indXYZT + t*indXYZ + k*indXY + j*np;
                    indT3 = r*indXYZ  +            k*indXY + j*np;
                    for (i=0; i<np; i++) {
                        indIMGS = indT2 + i;
                        indIMG  = indT1 + i;
                        indS    = indT3 + i;
                        pIMGSrd[indIMGS] = pIMGrd[indIMG]*pSrd[indS] - pIMGid[indIMG]*pSid[indS];
                        pIMGSid[indIMGS] = pIMGrd[indIMG]*pSid[indS] + pIMGid[indIMG]*pSrd[indS];
                    }
                }
                }
            }
        }
    }
    else {
        #pragma omp parallel for private(i,j,k,t,r,indT1,indT2,indT3,indIMG,indIMGS,indS)
        for (r=0; r<nr; r++) {
            for (t=0; t<nt; t++) {
                for (k=0; k<ns; k++) {
                for (j=0; j<nv; j++) {
                    indT1 =             t*indXYZ + k*indXY + j*np;
                    indT2 = r*indXYZT + t*indXYZ + k*indXY + j*np;
                    indT3 = r*indXYZ  +            k*indXY + j*np;
                    for (i=0; i<np; i++) {
                        indIMGS = indT2 + i;
                        indIMG  = indT1 + i;
                        indS    = indT3 + i;
                        pIMGSrf[indIMGS] = pIMGrf[indIMG]*pSrf[indS] - pIMGif[indIMG]*pSif[indS];
                        pIMGSif[indIMGS] = pIMGrf[indIMG]*pSif[indS] + pIMGif[indIMG]*pSrf[indS];
                    }
                }
                }
            }
        }
    }
    
    /* Free memory */
    mxFree(sizeOUT);
    
}
