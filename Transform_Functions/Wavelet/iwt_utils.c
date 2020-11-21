#ifdef TYPE1_DIM1
  #ifndef TYPE_REAL
    #define MAC_OPS \
      sum_xre += yare[ii] * ha[jj] + ydre[ii] * hd[jj];  \
      sum_xim += yaim[ii] * ha[jj] + ydim[ii] * hd[jj];  
  #else
    #define MAC_OPS \
      sum_xre += yare[ii] * ha[jj] + ydre[ii] * hd[jj];        
  #endif
#else
  #ifndef TYPE_REAL
    #define MAC_OPS \
      ii2 = ii*yincr;               \
      sum_xre += yare[ii2] * ha[jj] + ydre[ii2] * hd[jj];  \
      sum_xim += yaim[ii2] * ha[jj] + ydim[ii2] * hd[jj]; 
  #else
    #define MAC_OPS \
      ii2 = ii*yincr;               \
      sum_xre += yare[ii2] * ha[jj] + ydre[ii2] * hd[jj];        
  #endif 
#endif 


#ifdef TYPE1_DIM1
  #ifndef TYPE_REAL
    #define MAC_OPS_H \
      sum_xre += yare[jj] * ha[ii] + ydre[jj] * hd[ii];  \
      sum_xim += yaim[jj] * ha[ii] + ydim[jj] * hd[ii];  
  #else
    #define MAC_OPS_H \
      sum_xre += yare[jj] * ha[ii] + ydre[jj] * hd[ii];        
  #endif
#else
  #ifndef TYPE_REAL
    #define MAC_OPS_H \
      jj2 = jj*yincr;               \
      sum_xre += yare[jj2] * ha[ii] + ydre[jj2] * hd[ii];  \
      sum_xim += yaim[jj2] * ha[ii] + ydim[jj2] * hd[ii]; 
  #else
    #define MAC_OPS_H \
      jj2 = jj*yincr;               \
      sum_xre += yare[jj2] * ha[ii] + ydre[jj2] * hd[ii];        
  #endif 
#endif               

 
#ifdef TYPE_DIM1
  #define RECORD_OUTPUT_IDX kk2 = kk-crop2;
#else
  #define RECORD_OUTPUT_IDX kk2 = (kk-crop2)*xincr;
#endif   
  

#ifdef TYPE_REAL
  #define RECORD_OUTPUT \
    RECORD_OUTPUT_IDX;             \
    xre[kk2] = sum_xre;          
#else
  #define RECORD_OUTPUT \
    RECORD_OUTPUT_IDX;             \
    xre[kk2] = sum_xre;            \
    xim[kk2] = sum_xim;        
#endif  

                       
#ifdef TYPE_REAL
  #define INITIALIZE_SUMS sum_xre = 0.0; 
#else
  #define INITIALIZE_SUMS sum_xre = sum_xim = 0.0; 
#endif 



void UPCONV_FUNC_NAME(TYPE_NAME) 

#ifdef TYPE_DOUBLE
(double *yare, double *yaim, double *ydre, double *ydim, int ylen, int yincr, double *ha, double *hd, int hlen, double *xre, double *xim, int xincr, int crop2, int odd_dim)    
#else
(float *yare, float *yaim, float *ydre, float *ydim, int ylen, int yincr, float *ha, float *hd, int hlen, float *xre, float *xim, int xincr, int crop2, int odd_dim)    
#endif

{
  int nhm1 = hlen-1;
  int nym1 = 2*ylen-1;
  int ncm1 = 2*ylen + hlen - 1 - crop2 - odd_dim;
  int kk, ii, jj, iistart, iistop, pp, ii2, kk2, jj2;  
  
  #ifdef TYPE_DOUBLE  
    double sum_xre, sum_xim;        
  #else
    float sum_xre, sum_xim;        
  #endif    
         
    
  if (ylen >= hlen)
  {
    pp = 1; /* assumes hlen is even */  
      
    /* start up */
    for (kk=crop2; kk<nhm1; kk++)
    {
      INITIALIZE_SUMS;
                  
      for (ii=0, jj=kk; jj>=0; ii++, jj-=2) { MAC_OPS; }

      RECORD_OUTPUT;
    }
      
    /* middle */
    /* for (; kk<=nym1; kk++) */
    for (; kk<=ncm1; kk++)
    {
      INITIALIZE_SUMS;
                            
      iistart = (kk-nhm1+pp)/2;
      iistop = (kk-1+pp)/2;
      
      for (ii=iistart, jj=nhm1-pp; ii<=iistop; ii++, jj-=2) { MAC_OPS; }

      RECORD_OUTPUT;
      pp = 1-pp;
    }  
    return;
  
    /* end */
    for (; kk<=ncm1; kk++)
    {
      INITIALIZE_SUMS;
              
      iistart = (kk-nhm1+pp)/2;
      iistop = (kk-1+pp)/2;
    
      for (ii=iistart, jj=nhm1-pp; ii<=iistop; ii++, jj-=2) { MAC_OPS; }
    
      RECORD_OUTPUT;
      pp = 1-pp;
    }           
  }
  else /* ylen < hlen */ /* As of now, this section should never be reached since it doesn't properly reconstruct the signal */
  {
    /*  
    if ((ylen/2)*2==ylen) { pp = 1; }
    else                  { pp = 0; }           
     */
    
    pp = 1; /* ??? */   
      
    /* start up */
    for (kk=crop2; kk<nym1; kk++)
    {
      INITIALIZE_SUMS;
                  
      for (ii=0, jj=kk; jj>=0; ii++, jj-=2) { MAC_OPS_H; }

      RECORD_OUTPUT;
    }
      
    /* middle */
    /* for (; kk<=nhm1; kk++) */
    for (; kk<=ncm1; kk++)
    {
      INITIALIZE_SUMS
                            
      iistart = (kk-nym1+pp)/2;
      iistop = (kk-1+pp)/2;
      
      for (ii=iistart, jj=nym1-pp; ii<=iistop; ii++, jj-=2) { MAC_OPS_H }

      RECORD_OUTPUT;
      pp = 1-pp;
    }  
    return;
  
    /* end */
    for (; kk<=ncm1; kk++)
    {
      INITIALIZE_SUMS
              
      iistart = (kk-nym1+pp)/2;
      iistop = (kk-1+pp)/2;
    
      for (ii=iistart, jj=nym1-pp; ii<=iistop; ii++, jj-=2) { MAC_OPS_H }
    
      RECORD_OUTPUT;
      pp = 1-pp;
    }  
  }
}
    
 

          
#undef TYPE_NAME         
#undef MAC_OPS
#undef MAC_OPS_H
#undef RECORD_OUTPUT_IDX
#undef RECORD_OUTPUT     
#undef INITIALIZE_SUMS
          
          
          
     