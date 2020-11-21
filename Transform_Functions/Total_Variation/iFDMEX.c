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
  mxArray *I;
  mwSize  np, nv, ns, nt, ind, sh, sh_rt, sh_nt, sh_ns, jump, dir;
  long long rt;
  double  *pdir, *pIr, *pIi, *pFDr, *pFDi;
  float   *pdirf, *pIrf, *pIif, *pFDrf, *pFDif;
  
  /* Get sizes */
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
	if (precision == mxDOUBLE_CLASS) {
        if (mxIsComplex(right[0])) {
            I = mxCreateNumericArray(nD,sz,precision,mxCOMPLEX);
            pFDi = mxGetPi(right[0]);
            pIi = mxGetPi(I);
        }
        else {
            I = mxCreateNumericArray(nD,sz,precision,mxREAL);
        }
        pFDr = mxGetPr(right[0]);		
		pIr = mxGetPr(I);
	}
	else {
        if (mxIsComplex(right[0])) {
            I = mxCreateNumericArray(nD,sz,precision,mxCOMPLEX);
            pFDif = mxGetImagData(right[0]);
            pIif = mxGetImagData(I);
        }
        else {
            I = mxCreateNumericArray(nD,sz,precision,mxREAL);
        }
		pFDrf = mxGetData(right[0]);
		pIrf = mxGetData(I);
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
  
  /* Compute finite differences */
	if (precision == mxDOUBLE_CLASS) {
		
		if (dir == 1) {
			jump = -1;
			
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
                                for (np=1; np<NP; np++) {
                                    ind = sh + np;
                                    pIr[ind] = pFDr[ind+jump] - pFDr[ind];
                                    pIi[ind] = pFDi[ind+jump] - pFDi[ind];
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
                                for (np=1; np<NP; np++) {
                                    ind = sh + np;
                                    pIr[ind] = pFDr[ind+jump] - pFDr[ind];
                                }
                                
                                /*ind++;
                                 * pFDr[ind] = -pIr[ind];
                                 * pFDi[ind] = -pIi[ind];*/
                                
                            }
                        }
                    }
                }
            }
		}
		else if (dir == 2) {
			jump = -NP;
			
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=1; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIr[ind] = pFDr[ind+jump] - pFDr[ind];
                                    pIi[ind] = pFDi[ind+jump] - pFDi[ind];
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
            else {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=1; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIr[ind] = pFDr[ind+jump] - pFDr[ind];
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
		}
		
		else if (dir == 3) {
			jump = -NP*NV;
			
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=1; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIr[ind] = pFDr[ind+jump] - pFDr[ind];
                                    pIi[ind] = pFDi[ind+jump] - pFDi[ind];
                                }
                            }
                        }
                        
                        /*sh_ns = (NS-1)* NV*NP;
                         * for (nv=0; nv<NV; nv++) {
                         * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                         * for (np=0; np<NP; np++) {
                         * ind = sh + np;
                         * pIr[ind] = -pFDr[ind];
                         * pIi[ind] = -pFDi[ind];
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
                        for (ns=1; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIr[ind] = pFDr[ind+jump] - pFDr[ind];
                                }
                            }
                        }
                        
                        /*sh_ns = (NS-1)* NV*NP;
                         * for (nv=0; nv<NV; nv++) {
                         * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                         * for (np=0; np<NP; np++) {
                         * ind = sh + np;
                         * pIr[ind] = -pFDr[ind];
                         * pIi[ind] = -pFDi[ind];
                         * }
                         * }*/
                        
                    }
                }
            }
		}
		
		else if (dir == 4) {
			jump = -NP*NV*NS;
			
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=1; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIr[ind] = pFDr[ind+jump] - pFDr[ind];
                                    pIi[ind] = pFDi[ind+jump] - pFDi[ind];
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
                     * pIr[ind] = -pFDr[ind];
                     * pIi[ind] = -pFDi[ind];
                     * }
                     * }
                     * }*/
                    
                }
            }
            else {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=1; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIr[ind] = pFDr[ind+jump] - pFDr[ind];
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
                     * pIr[ind] = -pFDr[ind];
                     * pIi[ind] = -pFDi[ind];
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
			jump = -1;
			
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
                                for (np=1; np<NP; np++) {
                                    ind = sh + np;
                                    pIrf[ind] = pFDrf[ind+jump] - pFDrf[ind];
                                    pIif[ind] = pFDif[ind+jump] - pFDif[ind];
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
                                for (np=1; np<NP; np++) {
                                    ind = sh + np;
                                    pIrf[ind] = pFDrf[ind+jump] - pFDrf[ind];
                                }
                                
                                /*ind++;
                                 * pFDr[ind] = -pIr[ind];
                                 * pFDi[ind] = -pIi[ind];*/
                                
                            }
                        }
                    }
                }
            }
		}
		else if (dir == 2) {
			jump = -NP;
			
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=1; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIrf[ind] = pFDrf[ind+jump] - pFDrf[ind];
                                    pIif[ind] = pFDif[ind+jump] - pFDif[ind];
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
            else {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=1; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIrf[ind] = pFDrf[ind+jump] - pFDrf[ind];
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
		}
		
		else if (dir == 3) {
			jump = -NP*NV;
			
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=0; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=1; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIrf[ind] = pFDrf[ind+jump] - pFDrf[ind];
                                    pIif[ind] = pFDif[ind+jump] - pFDif[ind];
                                }
                            }
                        }
                        
                        /*sh_ns = (NS-1)* NV*NP;
                         * for (nv=0; nv<NV; nv++) {
                         * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                         * for (np=0; np<NP; np++) {
                         * ind = sh + np;
                         * pIr[ind] = -pFDr[ind];
                         * pIi[ind] = -pFDi[ind];
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
                        for (ns=1; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIrf[ind] = pFDrf[ind+jump] - pFDrf[ind];
                                }
                            }
                        }
                        
                        /*sh_ns = (NS-1)* NV*NP;
                         * for (nv=0; nv<NV; nv++) {
                         * sh = sh_rt + sh_nt + sh_ns + nv*NP;
                         * for (np=0; np<NP; np++) {
                         * ind = sh + np;
                         * pIr[ind] = -pFDr[ind];
                         * pIi[ind] = -pFDi[ind];
                         * }
                         * }*/
                        
                    }
                }
            }
		}
		
		else if (dir == 4) {
			jump = -NP*NV*NS;
			
            if (mxIsComplex(right[0])) {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=1; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIrf[ind] = pFDrf[ind+jump] - pFDrf[ind];
                                    pIif[ind] = pFDif[ind+jump] - pFDif[ind];
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
                     * pIr[ind] = -pFDr[ind];
                     * pIi[ind] = -pFDi[ind];
                     * }
                     * }
                     * }*/
                    
                }
            }
            else {
                #pragma omp parallel for private(rt, sh_rt, nt, sh_nt, ns, sh_ns, nv, sh, np, ind) shared(jump)
                for (rt=0; rt<RT; rt++) {
                    sh_rt = rt * NT*NS*NV*NP;
                    for (nt=1; nt<NT; nt++) {
                        sh_nt = nt * NS*NV*NP;
                        for (ns=0; ns<NS; ns++) {
                            sh_ns = ns * NV*NP;
                            for (nv=0; nv<NV; nv++) {
                                sh = sh_rt + sh_nt + sh_ns + nv*NP;
                                for (np=0; np<NP; np++) {
                                    ind = sh + np;
                                    pIrf[ind] = pFDrf[ind+jump] - pFDrf[ind];
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
                     * pIr[ind] = -pFDr[ind];
                     * pIi[ind] = -pFDi[ind];
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
  left[0] = I;

}
