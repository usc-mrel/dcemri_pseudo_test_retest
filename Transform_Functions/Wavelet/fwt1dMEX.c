/***********************
 * 1D forward wavelet transform mex file
 * - can operate on any dimension of N-dimensional data
 * - threaded with pthreads
 * 
 * - real or complex data, even-length filters, single- or double-precision
 *
 * [A,D] = fwt1dMEX(X,ha,hd,dim)
 * 
 * Travis Smith
 * October, 2011
 *
 ***********************/

#include "mex.h"
#include "math.h"
#include "string.h"
#include "pthread.h"



#include "wt_utils.h"
#include "fast_mxArray_setup.c"
#include "nd_indexing_setup.c"

#define _GNU_SOURCE_
#include "sched.h"

#ifndef MAXCORES
  #define MAXCORES 1
#endif       
        
/* turns on cpu affinity control */        
/*        
#define SETAFFINITY        
*/
        
/* prints some debug info */
/*
#define VERBOSE
*/         
             
        
typedef struct
{
  mwSize M_H, M_X, M_AD;
  mwSize *start_X, *start_AD, incr_X, incr_AD, nstart_X, nstart_AD;
    
  double *pXr_d, *pXi_d, *pAr_d, *pAi_d, *pDr_d, *pDi_d, *pHA_d, *pHD_d;
  float  *pXr_f, *pXi_f, *pAr_f, *pAi_f, *pDr_f, *pDi_f, *pHA_f, *pHD_f;
  
  mxClassID precision;
  mxComplexity cmplx;
  int extmode;
} array_params;
        
typedef struct 
{
  int id;
  mwSize idx_start, idx_stop;
  array_params *arrays;
} core_params;
    

/* interface to convolution functions */
void *firdn2(void *params)
{ 
  core_params *p = (core_params *)params; 
  array_params *a = p->arrays;
  int mm;  
    
  void (*fconv_d)(double *,double *,int,int,double *,double *,int,double *,double *,double *,double *,int,int) = NULL;
  void (*fconv_f)(float *,float *,int,int,float *,float *,int,float *,float *,float *,float *,int,int) = NULL;
    
  #ifdef SETAFFINITY
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(p->id*2,&set);
  sched_setaffinity(0,sizeof(set),&set);
  #endif   
             
  if (a->incr_X==1) { if (a->cmplx!=mxCOMPLEX) { if (a->precision==mxDOUBLE_CLASS) { fconv_d = &convdwn_1rd; }
                                                 else                              { fconv_f = &convdwn_1rf; }}
                      else                     { if (a->precision==mxDOUBLE_CLASS) { fconv_d = &convdwn_1cd; }
                                                 else                              { fconv_f = &convdwn_1cf; }}}
  else              { if (a->cmplx!=mxCOMPLEX) { if (a->precision==mxDOUBLE_CLASS) { fconv_d = &convdwn_2rd; }
                                                 else                              { fconv_f = &convdwn_2rf; }}
                      else                     { if (a->precision==mxDOUBLE_CLASS) { fconv_d = &convdwn_2cd; }
                                                 else                              { fconv_f = &convdwn_2cf; }}}                                  
                                     
  if (a->precision==mxDOUBLE_CLASS)
  {
    for (mm=p->idx_start; mm<=p->idx_stop; mm++)    
    {           
      fconv_d(a->pXr_d+a->start_X[mm],a->pXi_d+a->start_X[mm], a->M_X, a->incr_X, a->pHA_d, a->pHD_d, a->M_H, a->pAr_d+a->start_AD[mm], a->pAi_d+a->start_AD[mm], a->pDr_d+a->start_AD[mm], a->pDi_d+a->start_AD[mm], a->incr_AD, a->extmode);
    }
  }
  else
  {
    for (mm=p->idx_start; mm<=p->idx_stop; mm++)    
    {           
      fconv_f(a->pXr_f+a->start_X[mm],a->pXi_f+a->start_X[mm], a->M_X, a->incr_X, a->pHA_f, a->pHD_f, a->M_H, a->pAr_f+a->start_AD[mm], a->pAi_f+a->start_AD[mm], a->pDr_f+a->start_AD[mm], a->pDi_f+a->start_AD[mm], a->incr_AD, a->extmode);
    }  
  }                                 
                      
}

      
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{  
  /* Declare variables */
  const mwSize *size_X;
  mwSize *size_AD, ndims, opdim, iters_per_core, cc, NE=1; 
  pthread_t threadid[MAXCORES];
  array_params arrays;
  core_params params[MAXCORES];
  int ncores;
  pthread_attr_t attr;
  
  #ifdef SETAFFINITY
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(2*(NCORES-1),&set);
  sched_setaffinity(0,sizeof(set),&set);
  #endif
          
  pthread_attr_init(&attr);    
  
  /* access the input data */
  if (nrhs!=4) { mexErrMsgTxt("Need 4 input arguments."); } 
  
  arrays.precision = mxGetClassID(prhs[0]);
  
  if (arrays.precision==mxDOUBLE_CLASS)
  {
    arrays.pXr_d = (double *)mxGetData(prhs[0]);
    arrays.pXi_d = (double *)mxGetImagData(prhs[0]);
    arrays.pHA_d = (double *)mxGetData(prhs[1]);
    arrays.pHD_d = (double *)mxGetData(prhs[2]);  
    arrays.cmplx = ((arrays.pXi_d==NULL) ? mxREAL : mxCOMPLEX);
  }
  else
  {
    arrays.pXr_f = (float *)mxGetData(prhs[0]);
    arrays.pXi_f = (float *)mxGetImagData(prhs[0]);
    arrays.pHA_f = (float *)mxGetData(prhs[1]);
    arrays.pHD_f = (float *)mxGetData(prhs[2]); 
    arrays.cmplx = ((arrays.pXi_f==NULL) ? mxREAL : mxCOMPLEX);
  }
             
  opdim = *mxGetPr(prhs[3])-1;
  arrays.extmode = 2; /* 0=no extenstion, 1=constant extension, 2 = half-point symmetric extension */  
  
  /* get size of input data and set size of output data */
  arrays.M_H = mxGetM(prhs[1]);
  
  ndims = mxGetNumberOfDimensions(prhs[0]);
  if (opdim+1 > ndims) { mexErrMsgTxt("array does not have that many dimensions."); }
  size_X = mxGetDimensions(prhs[0]);
  for (cc=0; cc<ndims; cc++) { NE *= size_X[cc]; }
  
  arrays.M_X = size_X[opdim]; /* size of the dimension that we're filtering along */  
  arrays.M_AD = (arrays.M_X + arrays.M_H - 1)*0.5; 
  
  size_AD = (mwSize *)mxMalloc(ndims*sizeof(mwSize));
  memcpy(size_AD,size_X,ndims*sizeof(mwSize));
  size_AD[opdim] = arrays.M_AD; 

  /* setup outputs */
  if (arrays.precision==mxDOUBLE_CLASS)
  {
    create_array_d(&(plhs[0]), &arrays.pAr_d, &arrays.pAi_d, ndims, size_AD, arrays.cmplx, 0);
    create_array_d(&(plhs[1]), &arrays.pDr_d, &arrays.pDi_d, ndims, size_AD, arrays.cmplx, 0);
  }
  else
  {
    create_array_f(&(plhs[0]), &arrays.pAr_f, &arrays.pAi_f, ndims, size_AD, arrays.cmplx, 0);
    create_array_f(&(plhs[1]), &arrays.pDr_f, &arrays.pDi_f, ndims, size_AD, arrays.cmplx, 0);  
  }  

  /* setup n-dim indexing */
  nd_setup((mwSize *)size_X,ndims,opdim,&arrays.start_X,&arrays.nstart_X, &arrays.incr_X);
  nd_setup((mwSize *)size_AD,ndims,opdim,&arrays.start_AD,&arrays.nstart_AD, &arrays.incr_AD);
  
  /* setup param structs to pass to threads */
  /*
  ncores = MAXCORES;
  iters_per_core = arrays.nstart_X/ncores;
  if (iters_per_core==0 || (opdim==0 && NE<=8192)) { ncores = 1; }
  */ 
  
  ncores = ceil((float)NE / 8192.0);
  ncores = ((NE < 65536.0) ? 1 : ncores);
  ncores = ((ncores < 1) ? 1 : ncores);
  ncores = ((ncores > MAXCORES) ? MAXCORES : ncores);
  iters_per_core = ceil((float)arrays.nstart_X/(float)ncores);

  for (cc=0; cc<ncores; cc++)
  {
    params[cc].arrays = &arrays;
    params[cc].idx_start = cc*iters_per_core;
    params[cc].idx_stop = (cc+1)*iters_per_core-1;
    params[cc].id = cc;
  }
  params[ncores-1].idx_stop = arrays.nstart_X-1;

  #ifdef VERBOSE
  printf("MAXCORES = %d, ncores = %d, NE = %d, # iters in each thread: ",MAXCORES,ncores,NE);
  for (cc=0; cc<ncores; cc++) { printf("%d ",params[cc].idx_stop-params[cc].idx_start+1); }
  printf("\n");
  #endif  
  
  /* perform convolutions */
  for (cc=0; cc<ncores-1; cc++)
  {
    pthread_create(&threadid[cc],&attr,firdn2,(void *)(&params[cc]));
  }
  firdn2((void *)(&params[ncores-1]));
    
  pthread_attr_destroy(&attr);
  for (cc=0; cc<ncores-1; cc++)
  {
    pthread_join(threadid[cc], NULL); 
  }
  
  /* free memory */
  mxFree(arrays.start_X); 
  mxFree(arrays.start_AD);
  mxFree(size_AD);                  
  
}
