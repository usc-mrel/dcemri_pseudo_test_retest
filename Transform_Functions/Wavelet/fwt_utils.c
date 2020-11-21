
#ifdef TYPE1_DIM1
  #ifndef TYPE_REAL
    #define MAC_OPS \
      sum_are += xre[ii] * ha[jj];  \
      sum_aim += xim[ii] * ha[jj];  \
      sum_dre += xre[ii] * hd[jj];  \
      sum_dim += xim[ii] * hd[jj];  
  #else
    #define MAC_OPS \
      sum_are += xre[ii] * ha[jj];  \
      sum_dre += xre[ii] * hd[jj];         
  #endif
#else
  #ifndef TYPE_REAL
    #define MAC_OPS \
      ii2 = ii*xincr;               \
      sum_are += xre[ii2] * ha[jj];  \
      sum_aim += xim[ii2] * ha[jj];  \
      sum_dre += xre[ii2] * hd[jj];  \
      sum_dim += xim[ii2] * hd[jj];  
  #else
    #define MAC_OPS \
      ii2 = ii*xincr;               \
      sum_are += xre[ii2] * ha[jj];  \
      sum_dre += xre[ii2] * hd[jj];         
  #endif 
#endif              
              

#ifdef TYPE1_DIM1
  #ifndef TYPE_REAL
    #define MAC_OPS_H \
      sum_are += xre[jj] * ha[ii];  \
      sum_aim += xim[jj] * ha[ii];  \
      sum_dre += xre[jj] * hd[ii];  \
      sum_dim += xim[jj] * hd[ii];  
  #else
    #define MAC_OPS_H \
      sum_are += xre[jj] * ha[ii];  \
      sum_dre += xre[jj] * hd[ii];         
  #endif
#else
  #ifndef TYPE_REAL
    #define MAC_OPS_H \
      jj2 = jj*xincr;               \
      sum_are += xre[jj2] * ha[ii];  \
      sum_aim += xim[jj2] * ha[ii];  \
      sum_dre += xre[jj2] * hd[ii];  \
      sum_dim += xim[jj2] * hd[ii];  
  #else
    #define MAC_OPS_H \
      jj2 = jj*xincr;               \
      sum_are += xre[jj2] * ha[ii];  \
      sum_dre += xre[jj2] * hd[ii];         
  #endif 
#endif               
 
              
#ifdef TYPE_DIM1
  #define RECORD_OUTPUT_IDX kkm2 = kk/2;
#else
  #define RECORD_OUTPUT_IDX kkm2 = (kk/2)*yincr;
#endif  
                  
              
#ifdef TYPE_REAL
  #define RECORD_OUTPUT \
    RECORD_OUTPUT_IDX;             \
    yare[kkm2] = sum_are;         \
    ydre[kkm2] = sum_dre;           
#else
  #define RECORD_OUTPUT \
    RECORD_OUTPUT_IDX;             \
    yare[kkm2] = sum_are;         \
    yaim[kkm2] = sum_aim;         \
    ydre[kkm2] = sum_dre;         \
    ydim[kkm2] = sum_dim;  
#endif      
 


#ifdef TYPE_REAL
  #define INITIALIZE_SUMS sum_are = sum_dre = 0.0;
#else
  #define INITIALIZE_SUMS sum_are = sum_dre = sum_aim = sum_dim = 0.0;
#endif  
    
    
    

    

 
void CONVDWN_FUNC_NAME(TYPE_NAME) 

#ifdef TYPE_DOUBLE
(double *xre, double *xim, int xlen, int xincr, double *ha, double *hd, int hlen, double *yare, double *yaim, double *ydre, double *ydim, int yincr, int extmode)
#else
(float *xre, float *xim, int xlen, int xincr, float *ha, float *hd, int hlen, float *yare, float *yaim, float *ydre, float *ydim, int yincr, int extmode)
#endif

{
  int nhm1 = hlen-1;
  int nxm1 = xlen-1;
  int ncm1 = xlen + hlen - 2;
  int kk, ii, jj, kkm2, ii2, jj2;  
  int extincr=1;
  
  #ifdef TYPE_DOUBLE  
    double sum_are, sum_aim, sum_dre, sum_dim;      
  #else
    float sum_are, sum_aim, sum_dre, sum_dim;      
  #endif    
  
  if (extmode==1) { extincr = 0; }
 
  if (xlen >= hlen)
  {          
    /* start up */   
    for (kk=1; kk<=nhm1; kk+=2)  
    {      
      INITIALIZE_SUMS;          
      
      if (extmode>0) { for (ii=0, jj=kk+1; jj<=nhm1; jj++, ii+=extincr) { MAC_OPS; }} /* extension */                                                  
      for (ii=0, jj=kk; ii<=kk; ii++, jj--) { MAC_OPS; }      
      
      RECORD_OUTPUT;
    }     
        
    /* middle */
    for (; kk<=nxm1; kk+=2)
    {
      INITIALIZE_SUMS;
              
      for (ii=kk-nhm1, jj=nhm1; ii<=kk; ii++, jj--) { MAC_OPS; }
      
      RECORD_OUTPUT;
    }   
 
    /* end */
    for (; kk<=ncm1; kk+=2)
    {  
      INITIALIZE_SUMS;
              
      for (ii=kk-nhm1, jj=nhm1; ii<=nxm1; ii++, jj--) { MAC_OPS; }              
      if (extmode>0) { for (ii=nxm1; jj>=0; jj--, ii-=extincr) { MAC_OPS; }} /* extension */
    
      RECORD_OUTPUT;
    }  
                                          
  }
  else /* xlen < hlen */        
  {
          
    /* start up */   
    for (kk=1; kk<=nxm1; kk+=2)  
    {
      INITIALIZE_SUMS;
          
      if (extmode>0) { for (ii=0, jj=kk+1; jj<=nxm1; jj++, ii+=extincr) { MAC_OPS_H; }} /* extension */                                                  
      for (ii=0, jj=kk; ii<=kk; ii++, jj--) { MAC_OPS_H; }
      
      RECORD_OUTPUT;
    }     
        
    /* middle */
    for (; kk<=nhm1; kk+=2)
    {
      INITIALIZE_SUMS;
              
      for (ii=kk-nxm1, jj=nxm1; ii<=kk; ii++, jj--) { MAC_OPS_H; }
      
      RECORD_OUTPUT;
    }   
  
    /* end */
    for (; kk<=ncm1; kk+=2)
    {  
      INITIALIZE_SUMS;
              
      for (ii=kk-nxm1, jj=nxm1; ii<=nhm1; ii++, jj--) { MAC_OPS_H; }              
      if (extmode>0) { for (ii=nhm1; jj>=0; jj--, ii-=extincr) { MAC_OPS_H; }} /* extension */
    
      RECORD_OUTPUT;
    }                                    
  }          
}
          
#undef TYPE_NAME   
#undef MAC_OPS
#undef MAC_OPS_H
#undef RECORD_OUTPUT_IDX
#undef RECORD_OUTPUT
#undef INITIALIZE_SUMS
          
          
          
          
          
     