
/* create_array_d and create_array_f are alternatives to mxCreateNumericArray()
 * in that they create mxArrays (and pointers to them).  Most importantly, 
 * these functions allow fast mxArray allocation without initializing each element
 * to zero. 
 *
 * Travis Smith
 * Summer, 2011
 *
 *
 * create_array_d is for double-precision data
 * create_array_f is for single-precision data
 *
 * IN:
 * ndims = the desired # of dimensions for mxarr (and the length of sz)
 * sz = the array of length ndims which specifies the number of elements in each dimension
 * cx = mxREAL or mxCOMPLEX
 * z = 1 to zero each element in mxarr, = 0 to skip this time-consuming step (note: mxarr will be filled with garbage if z=0)
 * OUT:
 * mxarr = pointer to mxArray pointer (usually this is &(plhs[ii]))
 * pr, pi = pointers to {double,float} pointers (which can be used to access values in mxarr) 
 *
 * note: you can pass pi as NULL if cx=mxREAL
 *
 */

void create_array_d(mxArray **mxarr, double **pr, double **pi, mwSize ndims, mwSize *sz, mxComplexity cx, int z)
{
  int ii,numel,nbytes,sb = sizeof(double);
  
  *mxarr = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,cx);
  
  mxSetDimensions(*mxarr,sz,ndims);        
  numel = sz[0];        
  for (ii=1; ii<ndims; ii++) { numel *= sz[ii]; }        
      
  if (z==0)
  {
    nbytes = numel*sb;
    *pr = (double *)mxGetData(*mxarr);
    *pr = (double *)mxRealloc(*pr,nbytes);
    mxSetData(*mxarr,*pr);
    if (cx==mxCOMPLEX) 
    { 
      *pi = (double *)mxGetImagData(*mxarr);  
      *pi = (double *)mxRealloc(*pi,nbytes); 
      mxSetImagData(*mxarr,*pi);
    }
  }
  else
  {
    mxFree((void *)*pr);
    *pr = (double *)mxGetData(*mxarr);
    *pr = (double *)mxCalloc(numel,sb);
    mxSetData(*mxarr,*pr);
    if (cx==mxCOMPLEX) 
    { 
      *pi = (double *)mxGetImagData(*mxarr);  
      mxFree((void *)*pi);  
      *pi = (double *)mxCalloc(numel,sb); 
      mxSetImagData(*mxarr,*pi);
    }
  }      
}

void create_array_f(mxArray **mxarr, float **pr, float **pi, mwSize ndims, mwSize *sz, mxComplexity cx, int z)
{
  int ii,numel,nbytes,sb = sizeof(float);
  
  *mxarr = mxCreateNumericMatrix(0,0,mxSINGLE_CLASS,cx);  
         
  mxSetDimensions(*mxarr,sz,ndims);        
  numel = sz[0];        
  for (ii=1; ii<ndims; ii++) { numel *= sz[ii]; }        
      
  if (z==0)
  {
    nbytes = numel*sb;  
    *pr = (float *)mxGetData(*mxarr);
    *pr = (float *)mxRealloc(*pr,nbytes);
    mxSetData(*mxarr,*pr);
    if (cx==mxCOMPLEX) 
    { 
      *pi = (float *)mxGetImagData(*mxarr);  
      *pi = (float *)mxRealloc(*pi,nbytes); 
      mxSetImagData(*mxarr,*pi);
    }
  }
  else
  {
    *pr = (float *)mxGetData(*mxarr);  
    mxFree((void *)*pr);
    *pr = (float *)mxCalloc(numel,sb);
    mxSetData(*mxarr,*pr);
    if (cx==mxCOMPLEX) 
    { 
      *pi = (float *)mxGetImagData(*mxarr);  
      mxFree((void *)*pi);  
      *pi = (float *)mxCalloc(numel,sb); 
      mxSetImagData(*mxarr,*pi);
    }
  }    
}