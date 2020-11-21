/**************************************************************************
 MEX function to compute the approximate gradient of the absolute value
 
 Author: R. Marc Lebel
 Contact: mlebel@gmail.com
 Date: 12/2013
 
 Useage: k = convVS(img,sens)
 
 Input:
 img: numeric array (single/double; real/complex)
 sens: sensitivity maps
 
 Output:
 imgS: multichanel image
**************************************************************************/

#include "mex.h"
#include <math.h>
#include <string.h>
/*#include "fast_mxArray_setup.c"*/
/*#define FW 1.5*/ /* Standard deviation for gaussian weighting */
/*#define TR 4*/   /* Range of temporal convolution (+- this value) */

#ifdef __GNU__
    #include <omp.h>
#endif

#ifndef MAXCORES
  #define MAXCORES 1
#endif 


void mexFunction(int nlhs, mxArray *left[], int nrhs, const mxArray *right[]) {
    
    /*  Declare variables */
    int i, j, k, pe, r, t1;
    int TR;
    double FW;
    mwSize np, nv, ns, nt, nr, np_nv, np_ns, np_nt;
    mwSize nPE;
    mwSize indB, indBI, indI, indBO, indO;
    int pey, pez, pet, tmS, tmE;
    mwSize  nD, *sizeOUT;
	mxClassID precision;
    double  Gexpd, WTd;
    float   Gexpf, WTf;
    double  *pSZd, *pIMGrd, *pIMGid, *pUrd;
    double  *pIMGVSrd, *pIMGVSid, *pUVSrd;
	float   *pSZf, *pIMGrf, *pIMGif;
    float   *pIMGVSrf, *pIMGVSif, *pUVSrf;
    
    /*  Get dimensions, image size, and number of phase encodes */
    nD = 5;
    nPE = mxGetNumberOfElements(right[1])/3;
    if (mxGetClassID(right[2]) == mxDOUBLE_CLASS) {
        pSZd = mxGetPr(right[2]);
        np = (int)pSZd[0]; nv = (int)pSZd[1]; ns = (int)pSZd[2]; nt = (int)pSZd[3]; nr = (int)pSZd[4];
    }
    else {
        pSZf = mxGetData(right[2]);
        np = (int)pSZf[0]; nv = (int)pSZf[1]; ns = (int)pSZf[2]; nt = (int)pSZf[3]; nr = (int)pSZf[4];
    }
    /*mexPrintf("size: %i x %i x %i x %i x %i\n",np,nv,ns,nt,nr);
    mexPrintf("nPE: %i\n",nPE);*/
    
    /*  Get convolution size */
    TR = (int)mxGetScalar(right[3]);
    FW = mxGetScalar(right[4]);
    
    /*  Perform strange memory copy to replicate the size (needed for create_array_d/f) */
    sizeOUT = (mwSize *)mxMalloc(nD*sizeof(mwSize));
    sizeOUT[0] = np; sizeOUT[1] = nv; sizeOUT[2] = ns; sizeOUT[3] = nt; sizeOUT[4] = nr;
    /*mexPrintf("sizeOUT: %i x %i x %i x %i x %i\n",sizeOUT[0],sizeOUT[1],sizeOUT[2],sizeOUT[3],sizeOUT[4]);*/
    
    /*  Test for complex and obtain data class */
    if (!(mxIsComplex(right[0])))
        mexErrMsgTxt("Input image needs to be complex");
    if (mxIsComplex(right[1]))
        mexErrMsgTxt("Phase encode table cannot be complex");
    if (mxGetClassID(right[1]) == mxSINGLE_CLASS)
        mexErrMsgTxt("Phase encode table must be double");
   	precision = mxGetClassID(right[0]);
    
    
    /*  Create output arrays */
    left[0] = mxCreateNumericArray(nD,sizeOUT,precision,mxCOMPLEX);
    left[1] = mxCreateNumericArray(nD,sizeOUT,precision,mxREAL);
    
    /*  Get input/output pointers */
    pUrd   = mxGetPr(right[1]);
	if (precision == mxDOUBLE_CLASS) {
        /*  Input */
		pIMGrd = mxGetPr(right[0]);
        pIMGid = mxGetPi(right[0]);
        
		/*  Output */
        pIMGVSrd = mxGetPr(left[0]);
        pIMGVSid = mxGetPi(left[0]);
        pUVSrd   = mxGetPr(left[1]);
	}
	else {
        /*  Input */
		pIMGrf = mxGetData(right[0]);
        pIMGif = mxGetImagData(right[0]);
        
        /*  Output */
        pIMGVSrf = mxGetData(left[0]);
        pIMGVSif = mxGetImagData(left[0]);
        pUVSrf   = mxGetData(left[1]);
	}
    
    /*  Define common indices */
    np_nv = np*nv;
    np_ns = np*nv*ns;
    np_nt = np*nv*ns*nt;
    
    #ifdef __GNU__
        /* Set number of threads */
        omp_set_num_threads(MAXCORES);
	#endif
    
    /*  Loop through elements and perform convolution */
    if (precision == mxDOUBLE_CLASS) {
        Gexpd = -1.0/(2.0*FW*FW);
        #pragma omp parallel for private(i,j,k,pe,t1,r,pey,pez,pet,tmS,tmE,indB,indBI,indI,indBO,indO,WTd)
        for (pe=0; pe<nPE; pe++) {
            /*  Get current ky/kz/t location */
            pey = pUrd[pe]-1;
            pez = pUrd[nPE+pe]-1;
            pet = pUrd[2*nPE+pe]-1;
            
            /*  Define time convolution range */
            tmS = ((pet - TR) < 0) ? 0 : pet-TR;
            tmE = ((pet + TR) >= nt) ? nt : pet+TR;
            
            /*  Convolve this phase encode through time */
            for (t1=tmS; t1<tmE; t1++) {
                WTd = exp(Gexpd*(t1-pet)*(t1-pet))+0.01;
                
                for (i=0; i<np; i++) {
                    /*  Convert to base kx/ky/kz index */
                    indB = i + pey*np + pez*np_nv;
                    
                    /*  Define base input index */
                    indBI = indB + pet*np_ns;
                    indBO = i + pey*np + pez*np_nv + t1*np_ns;

                    for (r=0; r<nr; r++) {
                        
                        /*  Define input and output indices */
                        indI = indBI + r*np_nt;
                        indO = indBO + r*np_nt;
                        
                        /*  Perform actual convolution */
                        pIMGVSrd[indO] += pIMGrd[indI] * WTd;
                        pIMGVSid[indO] += pIMGid[indI] * WTd;
                        pUVSrd[indO] += WTd;
                    }
                }
            }
        }
        
        
        /*  Perform density correction */
        #pragma omp parallel for private(i,j,k,t1,r,indB,indBO,indO)
        for (i=0; i<np; i++) {
            for (j=0; j<nv; j++) {
            for (k=0; k<ns; k++) {
                indB = i + j*np + k*np_nv;
                for (t1=0; t1<nt; t1++){
                    indBO = indB + t1*np_ns;
                    for (r=0; r<nr; r++) {
                        indO = indBO + r*np_nt;
                        pIMGVSrd[indO] /= (pUVSrd[indO] + 1.0e-20);
                        pIMGVSid[indO] /= (pUVSrd[indO] + 1.0e-20);
                    }
                }
            }
            }
        }

    }
    else {
        Gexpf = -1.0/(2.0*FW*FW);
        #pragma omp parallel for private(i,j,k,pe,t1,r,pey,pez,pet,tmS,tmE,indB,indBI,indI,indBO,indO,WTf)
        for (pe=0; pe<nPE; pe++) {
            
            /*  Get current ky/kz/t location */
            pey = pUrd[pe]-1;
            pez = pUrd[nPE+pe]-1;
            pet = pUrd[2*nPE+pe]-1;
            
            /*  Define time convolution range */
            tmS = ((pet - TR) < 0) ? 0 : pet-TR;
            tmE = ((pet + TR) >= nt) ? nt : pet+TR;
            
            /*  Convolve this phase encode through time */
            for (t1=tmS; t1<tmE; t1++) {
                WTf = exp(Gexpf*(t1-pet)*(t1-pet)) + 0.01;
                
                for (i=0; i<np; i++) {
                    
                    /*  Convert to base kx/ky/kz index */
                    indB = i + pey*np + pez*np_nv;
                    
                    /*  Define base input and output indices */
                    indBI = indB + pet*np_ns;
                    indBO = i + pey*np + pez*np_nv + t1*np_ns;
                    
                    for (r=0; r<nr; r++) {
                        
                        /*  Define input and output indices */
                        indI = indBI + r*np_nt;
                        indO = indBO + r*np_nt;
                        
                        /*  Perform actual convolution */
                        pIMGVSrf[indO] += pIMGrf[indI] * WTf;
                        pIMGVSif[indO] += pIMGif[indI] * WTf;
                        pUVSrf[indO] += WTf;
                    }
                }
            }
        }
        
        
        /*  Perform density correction */
        #pragma omp parallel for private(i,j,k,t1,r,indB,indBO,indO)
        for (i=0; i<np; i++) {
            for (j=0; j<nv; j++) {
            for (k=0; k<ns; k++) {
                indB = i + j*np + k*np_nv;
                for (t1=0; t1<nt; t1++){
                    indBO = indB + t1*np_ns;
                    for (r=0; r<nr; r++) {
                        indO = indBO + r*np_nt;
                        pIMGVSrf[indO] /= (pUVSrf[indO] + 1.0e-20);
                        pIMGVSif[indO] /= (pUVSrf[indO] + 1.0e-20);
                    }
                }
            }
            }
        }
    }
    
    
    /* Free memory */
    mxFree(sizeOUT);
    
}
