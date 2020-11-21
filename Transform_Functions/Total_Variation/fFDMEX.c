#include "mex.h"

#ifdef __GNU__
    #include <omp.h>
#endif

#ifndef MAXCORES
  #define MAXCORES 1
#endif 

void mexFunction(int nlhs, mxArray *left[], int nrhs, const mxArray *right[]) {
    
    /* Declare variables */
    mwSize  nD, NP=1, NV=1, NS=1, NT=1, RT=1;
    const   mwSize *sz;
    mxClassID precision;
    mxArray *FD;
    mwSize  np, nv, ns, nt, ind, sh, sh_rt, sh_nt, sh_ns, jump, dir;
    long long rt;
    double  *pdir, *pIr, *pIi, *pFDr, *pFDi;
    float   *pdirf, *pIrf, *pIif, *pFDrf, *pFDif;
    
    /* Get sizes and data type */
    np = mxGetM(right[0]);
    nv = mxGetN(right[0]);
    nD = mxGetNumberOfDimensions(right[0]);
    sz = mxGetDimensions(right[0]);
    precision = mxGetClassID(right[0]);
    
    /* Determine how many rows, columns, depth, time, and other dimensions */
    /* This is done to parallelize the "other dimensions" */
    NP = sz[0];
    if (nD > 1)
        NV = sz[1];
    if (nD > 2)
        NS = sz[2];
    if (nD > 3)
        NT = sz[3];
    if (nD > 4) {
        for (np=4; np<nD; np++){
            RT *= sz[np];
        }
    }
    /*mexPrintf("NP: %i\n",NP);
    mexPrintf("NV: %i\n",NV);
    mexPrintf("NS: %i\n",NS);
    mexPrintf("NT: %i\n",NT);
    mexPrintf("RT: %i\n",RT);*/
    
    /*  Check if complex */
    /*if (!mxIsComplex(right[0]))
        mexErrMsgTxt("Function requires complex data");*/
    
    /* Create output and get input/output pointers */
    if (precision ==  mxDOUBLE_CLASS) {
        if (mxIsComplex(right[0])) {
            FD = mxCreateNumericArray(nD,sz,precision,mxCOMPLEX);
            pIi = mxGetPi(right[0]);
            pFDi = mxGetPi(FD);
        }
        else {
            FD = mxCreateNumericArray(nD,sz,precision,mxREAL);
        }
        pIr = mxGetPr(right[0]);
        pFDr = mxGetPr(FD);
    }
    else {
        if (mxIsComplex(right[0])) {
            FD = mxCreateNumericArray(nD,sz,precision,mxCOMPLEX);
            pIif = mxGetImagData(right[0]);
            pFDif = mxGetImagData(FD);
        }
        else {
            FD = mxCreateNumericArray(nD,sz,precision,mxREAL);
        }
        pIrf = mxGetData(right[0]);
        pFDrf = mxGetData(FD);
    }
    
    /*  Get transform direction */
    if (mxGetClassID(right[1]) == mxDOUBLE_CLASS) {
        pdir = mxGetData(right[1]);
        dir = (mwSize) pdir[0];
    }
    else {
        pdirf = mxGetData(right[1]);
        dir = (mwSize) pdirf[0];
    }
    
    #ifdef __GNU__
        /* Set number of threads */
        omp_set_num_threads(MAXCORES);
	#endif
    
    
    /*  Compute finite differences */
    if (precision == mxDOUBLE_CLASS) {
        
        if (dir == 1) {
            jump = 1;
            
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP-1; np++) {
                                    ind = sh + np;
                                    pFDr[ind] = pIr[ind+jump] - pIr[ind];
                                    pFDi[ind] = pIi[ind+jump] - pIi[ind];
                                }
                            
                                /*ind++;
                                * pFDr[ind] = -pIr[ind];
                                * pFDi[ind] = -pIi[ind];*/
                            
                            }
                        }
                    }
                }
            }
            else {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP-1; np++) {
                                    ind = sh + np;
                                    pFDr[ind] = pIr[ind+jump] - pIr[ind];
                                }
                            
                                /*ind++;
                                * pFDr[ind] = -pIr[ind];*/                            
                            }
                        }
                    }
                }
            }
        }
        
        
        else if (dir == 2) {
            jump = NP;
            
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV-1; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDr[ind] = pIr[ind+jump] - pIr[ind];
                                    pFDi[ind] = pIi[ind+jump] - pIi[ind];
                                }
                            }
                        
                            /*sh = sh_rt + sh_nt + sh_ns + (NV-1)*NP;
                            * for (np=0; np<NP; np++) {
                            * ind = sh + np;
                            * pFDr[ind] = -pIr[ind];
                            * pFDi[ind] = -pIi[ind];
                            * }*/
                        
                        }
                    }
                }
            }
            else{
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV-1; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDr[ind] = pIr[ind+jump] - pIr[ind];
                                }
                            }
                        
                            /*sh = sh_rt + sh_nt + sh_ns + (NV-1)*NP;
                            * for (np=0; np<NP; np++) {
                            * ind = sh + np;
                            * pFDr[ind] = -pIr[ind];
                            * }*/
                        
                        }
                    }
                }
            }
        }
        
        
        else if (dir == 3) {
            jump = NP*NV;
            
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS-1; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDr[ind] = pIr[ind+jump] - pIr[ind];
                                    pFDi[ind] = pIi[ind+jump] - pIi[ind];
                                }
                            }
                        }
                    
                        /*sh_ns = (NS-1)* NV*NP;
                        * for (nv=0; nv<NV; nv++) {
                        * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                        * for (np=0; np<NP; np++) {
                        * ind = sh + np;
                        * pFDr[ind] = -pIr[ind];
                        * }
                        * }*/
                    
                    }
                }
            }
            else {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS-1; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDr[ind] = pIr[ind+jump] - pIr[ind];
                                }
                            }
                        }
                    
                        /*sh_ns = (NS-1)* NV*NP;
                        * for (nv=0; nv<NV; nv++) {
                        * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                        * for (np=0; np<NP; np++) {
                        * ind = sh + np;
                        * pFDr[ind] = -pIr[ind];
                        * }
                        * }*/
                    
                    }
                }
            }
        }
        
        
        else if (dir == 4) {
            jump = NP*NV*NS;
            
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT-1; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDr[ind] = pIr[ind+jump] - pIr[ind];
                                    pFDi[ind] = pIi[ind+jump] - pIi[ind];
                                }
                            }
                        }
                    }
                
                    /*sh_nt = (NT-1)*NS*NV*NP;
                    * for (ns=0; ns<NS; ns++) {
                    * sh_ns = ns * NV*NP;
                    * for (nv=0; nv<NV; nv++) {
                    * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                    * for (np=0; np<NP; np++) {
                    * ind = sh + np;
                    * pFDr[ind] = -pIr[ind];
                    * }
                    * }
                    * }*/
                    
                }
            }
            else{
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT-1; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDr[ind] = pIr[ind+jump] - pIr[ind];
                                }
                            }
                        }
                    }
                
                    /*sh_nt = (NT-1)*NS*NV*NP;
                    * for (ns=0; ns<NS; ns++) {
                    * sh_ns = ns * NV*NP;
                    * for (nv=0; nv<NV; nv++) {
                    * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                    * for (np=0; np<NP; np++) {
                    * ind = sh + np;
                    * pFDr[ind] = -pIr[ind];
                    * }
                    * }
                    * }*/
                    
                }
            }
        }
        
        else
            mexErrMsgTxt("Unsupported transform direction");
        
    }
    
    /*  Single precision */
    else {
        
        if (dir == 1) {
            jump = 1;
            
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP-1; np++) {
                                    ind = sh + np;
                                    pFDrf[ind] = pIrf[ind+jump] - pIrf[ind];
                                    pFDif[ind] = pIif[ind+jump] - pIif[ind];
                                }
                            
                                /*ind++;
                                * pFDrf[ind] = -pIrf[ind];
                                */
                            
                            }
                        }
                    }
                }
            }
            else{
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP-1; np++) {
                                    ind = sh + np;
                                    pFDrf[ind] = pIrf[ind+jump] - pIrf[ind];
                                }
                            
                                /*ind++;
                                * pFDrf[ind] = -pIrf[ind];
                                */
                            
                            }
                        }
                    }
                }
            }
        }
        
        else if (dir == 2) {
            jump = NP;
            
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV-1; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDrf[ind] = pIrf[ind+jump] - pIrf[ind];
                                    pFDif[ind] = pIif[ind+jump] - pIif[ind];
                                }
                            }
                            
                            /*sh = sh_rt + sh_nt + sh_ns + (NV-1)*NP;
                            * for (np=0; np<NP; np++) {
                            * ind = sh + np;
                            * pFDrf[ind] = -pIrf[ind];
                            * pFDif[ind] = -pIif[ind];
                            * }*/
                        
                        }
                    }
                }
            }
            else{
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV-1; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDrf[ind] = pIrf[ind+jump] - pIrf[ind];
                                }
                            }
                            
                            /*sh = sh_rt + sh_nt + sh_ns + (NV-1)*NP;
                            * for (np=0; np<NP; np++) {
                            * ind = sh + np;
                            * pFDrf[ind] = -pIrf[ind];
                            * }*/
                        
                        }
                    }
                }
            }
        }
        
        else if (dir == 3) {
            jump = NP*NV;
            
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS-1; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDrf[ind] = pIrf[ind+jump] - pIrf[ind];
                                    pFDif[ind] = pIif[ind+jump] - pIif[ind];
                                }
                            }
                        }
                    
                        /*sh_ns = (NS-1)* NV*NP;
                        * for (nv=0; nv<NV; nv++) {
                        * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                        * for (np=0; np<NP; np++) {
                        * ind = sh + np;
                        * pFDrf[ind] = -pIrf[ind];
                        * pFDif[ind] = -pIif[ind];
                        * }
                        * }*/
                    
                    }
                }
            }
            else {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS-1; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDrf[ind] = pIrf[ind+jump] - pIrf[ind];
                                }
                            }
                        }
                    
                        /*sh_ns = (NS-1)* NV*NP;
                        * for (nv=0; nv<NV; nv++) {
                        * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                        * for (np=0; np<NP; np++) {
                        * ind = sh + np;
                        * pFDrf[ind] = -pIrf[ind];
                        * }
                        * }*/
                    
                    }
                }
            }
        }
        
        else if (dir == 4) {
            jump = NP*NV*NS;
            
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT-1; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDrf[ind] = pIrf[ind+jump] - pIrf[ind];
                                    pFDif[ind] = pIif[ind+jump] - pIif[ind];
                                }
                            }
                        }
                    }
                    
                    /*sh_nt = (NT-1)*NS*NV*NP;
                     * for (ns=0; ns<NS; ns++) {
                     * sh_ns = ns * NV*NP;
                     * for (nv=0; nv<NV; nv++) {
                     * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                     * for (np=0; np<NP; np++) {
                     * ind = sh + np;
                     * pFDrf[ind] = -pIrf[ind];
                     * pFDif[ind] = -pIif[ind];
                     * }
                     * }
                     * }*/
                
                }
            }
            else {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT-1; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pFDrf[ind] = pIrf[ind+jump] - pIrf[ind];
                                }
                            }
                        }
                    }
                    
                    /*sh_nt = (NT-1)*NS*NV*NP;
                     * for (ns=0; ns<NS; ns++) {
                     * sh_ns = ns * NV*NP;
                     * for (nv=0; nv<NV; nv++) {
                     * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                     * for (np=0; np<NP; np++) {
                     * ind = sh + np;
                     * pFDrf[ind] = -pIrf[ind];
                     * }
                     * }
                     * }*/
                
                }
            }
        }
        
        else
            mexErrMsgTxt("Unsupported transform direction");
        
    }
    
    /* Output */
    left[0] = FD;
    
}
