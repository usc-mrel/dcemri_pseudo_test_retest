

void create_array_d(mxArray **mxarr, double **pr, double **pi, mwSize ndims, mwSize *sz, mxComplexity cx, int z)
{
  int ii,numel,nbytes,sb = sizeof(double);
  
  *mxarr = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,cx);
  *pr = (double *)mxGetData(*mxarr);
  *pi = (double *)mxGetImagData(*mxarr);
  
  mxSetDimensions(*mxarr,sz,ndims);        
  numel = sz[0];        
  for (ii=1; ii<ndims; ii++) { numel *= sz[ii]; }        
      
  if (z==0)
  {
    nbytes = numel*sb;        
    *pr = (double *)mxRealloc(*pr,nbytes);
    mxSetData(*mxarr,*pr);
    if (cx==mxCOMPLEX) 
    { 
      *pi = (double *)mxRealloc(*pi,nbytes); 
      mxSetImagData(*mxarr,*pi);
    }
  }
  else
  {
    mxFree((void *)*pr);
    *pr = (double *)mxCalloc(numel,sb);
    mxSetData(*mxarr,*pr);
    if (cx==mxCOMPLEX) 
    { 
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
  *pr = (float *)mxGetData(*mxarr);
  *pi = (float *)mxGetImagData(*mxarr);
  
  mxSetDimensions(*mxarr,sz,ndims);        
  numel = sz[0];        
  for (ii=1; ii<ndims; ii++) { numel *= sz[ii]; }        
      
  if (z==0)
  {
    nbytes = numel*sb;        
    *pr = (float *)mxRealloc(*pr,nbytes);
    mxSetData(*mxarr,*pr);
    if (cx==mxCOMPLEX) 
    { 
      *pi = (float *)mxRealloc(*pi,nbytes); 
      mxSetImagData(*mxarr,*pi);
    }
  }
  else
  {
    mxFree((void *)*pr);
    *pr = (float *)mxCalloc(numel,sb);
    mxSetData(*mxarr,*pr);
    if (cx==mxCOMPLEX) 
    { 
      mxFree((void *)*pi);  
      *pi = (float *)mxCalloc(numel,sb); 
      mxSetImagData(*mxarr,*pi);
    }
  }    
}