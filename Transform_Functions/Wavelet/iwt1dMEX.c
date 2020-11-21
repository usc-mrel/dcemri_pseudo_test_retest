/***********************
 * 1D inverse wavelet transform mex file
 * - can operate on any dimension of N-dimensional data
 * - threaded with pthreads
 * 
 * - real or complex data, even-length filters, single- or double-precision
 *
 * [X] = iwt1dMEX(A,D,hlo,hhi,dim)
 * [X] = iwt1dMEX(A,D,hlo,hhi,dim,odd_dim)
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
  mwSize M_H, M_X, M_AD, crop, odd_dim;
  mwSize *start_X, *start_AD, incr_X, incr_AD, nstart_X, nstart_AD;
  
  double *pXr_d, *pXi_d, *pAr_d, *pAi_d, *pDr_d, *pDi_d, *pHA_d, *pHD_d;
  float  *pXr_f, *pXi_f, *pAr_f, *pAi_f, *pDr_f, *pDi_f, *pHA_f, *pHD_f;
  
  mxClassID cmplx, precision;
} array_params;
        
typedef struct 
{
  int id;
  mwSize idx_start, idx_stop;
  array_params *arrays;
} core_params;
     


/* interface to convolution functions */
void *upfir2(void *params)
{ 
  core_params *p = (core_params *)params; 
  array_params *a = p->arrays;
  int mm;
  int crop2 = ceil((float)a->crop/2.0);  
  
  void (*fconv_d)(double *,double *,double *,double *,int,int,double *,double *,int,double *,double *,int,int,int) = NULL; 
  void (*fconv_f)(float *,float *,float *,float *,int,int,float *,float *,int,float *,float *,int,int,int) = NULL;
  
  #ifdef SETAFFINITY
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(p->id*2,&set);
  sched_setaffinity(0,sizeof(set),&set);
  #endif   
  
  if (a->incr_X==1) { if (a->cmplx!=mxCOMPLEX) { if (a->precision==mxDOUBLE_CLASS) { fconv_d = &upconv_1rd; }
                                                 else                              { fconv_f = &upconv_1rf; }}
                      else                     { if (a->precision==mxDOUBLE_CLASS) { fconv_d = &upconv_1cd; }
                                                 else                              { fconv_f = &upconv_1cf; }}}
  else              { if (a->cmplx!=mxCOMPLEX) { if (a->precision==mxDOUBLE_CLASS) { fconv_d = &upconv_2rd; }
                                                 else                              { fconv_f = &upconv_2rf; }}
                      else                     { if (a->precision==mxDOUBLE_CLASS) { fconv_d = &upconv_2cd; }
                                                 else                              { fconv_f = &upconv_2cf; }}}           
          
  if (a->precision==mxDOUBLE_CLASS)
  {
    for (mm=p->idx_start; mm<=p->idx_stop; mm++)    
    {           
      fconv_d(a->pAr_d+a->start_AD[mm],a->pAi_d+a->start_AD[mm],a->pDr_d+a->start_AD[mm],a->pDi_d+a->start_AD[mm],a->M_AD,a->incr_AD,a->pHA_d,a->pHD_d,a->M_H,a->pXr_d+a->start_X[mm],a->pXi_d+a->start_X[mm],a->incr_X,crop2,a->odd_dim);
    }
  }
  else
  {
    for (mm=p->idx_start; mm<=p->idx_stop; mm++)    
    {           
      fconv_f(a->pAr_f+a->start_AD[mm],a->pAi_f+a->start_AD[mm],a->pDr_f+a->start_AD[mm],a->pDi_f+a->start_AD[mm],a->M_AD,a->incr_AD,a->pHA_f,a->pHD_f,a->M_H,a->pXr_f+a->start_X[mm],a->pXi_f+a->start_X[mm],a->incr_X,crop2,a->odd_dim);
    }  
  }      
}

        
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{  
  /* Declare variables */
  const mwSize *size_AD;
  mwSize *size_X, cc, ndims, opdim, iters_per_core, NE=1; 
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
  if (nrhs!=6 && nrhs!=5) { mexErrMsgTxt("Need 5 or 6 input arguments."); } 
  
  arrays.precision = mxGetClassID(prhs[0]);
  if (arrays.precision==mxDOUBLE_CLASS)
  {
    arrays.pAr_d = (double *)mxGetData(prhs[0]);
    arrays.pAi_d = (double *)mxGetImagData(prhs[0]);
    arrays.pDr_d = (double *)mxGetData(prhs[1]);
    arrays.pDi_d = (double *)mxGetImagData(prhs[1]);   
    arrays.pHA_d = (double *)mxGetData(prhs[2]);
    arrays.pHD_d = (double *)mxGetData(prhs[3]);     
    arrays.cmplx = ((arrays.pAi_d==NULL) ? mxREAL : mxCOMPLEX);   
  }
  else
  {
    arrays.pAr_f = (float *)mxGetData(prhs[0]);
    arrays.pAi_f = (float *)mxGetImagData(prhs[0]);
    arrays.pDr_f = (float *)mxGetData(prhs[1]);
    arrays.pDi_f = (float *)mxGetImagData(prhs[1]); 
    arrays.pHA_f = (float *)mxGetData(prhs[2]);
    arrays.pHD_f = (float *)mxGetData(prhs[3]);     
    arrays.cmplx = ((arrays.pAi_f==NULL) ? mxREAL : mxCOMPLEX);   
  }  
  
  opdim = *mxGetPr(prhs[4])-1;
  
  arrays.odd_dim = 0; /* if output will be odd in length, odd_dim should be 1 */
  if (nrhs==6) { arrays.odd_dim = *mxGetPr(prhs[5]); }
  
  /* get size of input data and set size of output data */
  arrays.M_H = mxGetM(prhs[2]);  
  
  ndims = mxGetNumberOfDimensions(prhs[0]);
  if (opdim+1 > ndims) { mexErrMsgTxt("array does not have that many dimensions."); }
  size_AD = mxGetDimensions(prhs[0]);
  for (cc=0; cc<ndims; cc++) { NE *= size_AD[cc]; }
  
  arrays.M_AD = size_AD[opdim]; /* size of dimension we're filtering in */
    
  arrays.crop = 2*arrays.M_H-3;    
  
  arrays.M_X = arrays.M_AD*2 + arrays.M_H - 1 - arrays.crop - arrays.odd_dim;
  
  if (arrays.M_H > arrays.M_X) { mexErrMsgTxt("cannot reconstruct: filter is shorter than original signal."); }
  
  size_X = (mwSize *)mxMalloc(ndims*sizeof(mwSize));
  memcpy(size_X,size_AD,ndims*sizeof(mwSize));
  size_X[opdim] = arrays.M_X; 
  
  /* setup outputs */ 
  if (arrays.precision==mxDOUBLE_CLASS)
  {
    create_array_d(&(plhs[0]), &arrays.pXr_d, &arrays.pXi_d, ndims, size_X, arrays.cmplx, 0);  
  }
  else
  {
    create_array_f(&(plhs[0]), &arrays.pXr_f, &arrays.pXi_f, ndims, size_X, arrays.cmplx, 0);
  }  
  
  /* setup n-dim indexing */
  nd_setup((mwSize *)size_X,ndims,opdim,&arrays.start_X,&arrays.nstart_X, &arrays.incr_X);
  nd_setup((mwSize *)size_AD,ndims,opdim,&arrays.start_AD,&arrays.nstart_AD, &arrays.incr_AD);
  
  /* setup param structs to pass to threads */
  /*
  ncores = MAXCORES;
  iters_per_core = arrays.nstart_AD/ncores;
  if (iters_per_core==0 || (opdim==0 && NE<=8192)) { ncores = 1; }
  */ 

  ncores = ceil((float)NE / 8192.0);
  ncores = ((NE < 65536.0) ? 1 : ncores);
  ncores = ((ncores < 1) ? 1 : ncores);
  ncores = ((ncores > MAXCORES) ? MAXCORES : ncores);
  iters_per_core = ceil((float)arrays.nstart_AD/(float)ncores);

  for (cc=0; cc<ncores; cc++)
  {
    params[cc].arrays = &arrays;
    params[cc].idx_start = cc*iters_per_core;
    params[cc].idx_stop = (cc+1)*iters_per_core-1;
    params[cc].id = cc;
  }
  params[ncores-1].idx_stop = arrays.nstart_AD-1;

#ifdef VERBOSE
  printf("MAXCORES = %d, ncores = %d, NE = %d, # iters in each thread: ",MAXCORES,ncores,NE);
  for (cc=0; cc<ncores; cc++) { printf("%d ",params[cc].idx_stop-params[cc].idx_start+1); }
  printf("\n");
#endif  
  
  /* perform convolutions */
  for (cc=0; cc<ncores-1; cc++)
  {
    pthread_create(&threadid[cc],&attr,upfir2,(void *)(&params[cc]));
  }
  upfir2((void *)(&params[ncores-1]));
    
  pthread_attr_destroy(&attr);
  for (cc=0; cc<ncores-1; cc++)
  {
    pthread_join(threadid[cc], NULL); 
  }
  
  /* free memory */
  mxFree(arrays.start_X); 
  mxFree(arrays.start_AD);
  mxFree(size_X);                  
  
}
