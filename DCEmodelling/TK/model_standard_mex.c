/*==========================================================
 * modelStandard.c - STARDCE toolbox
 *
 * Implements the DCE standard model
 *
 * The calling syntax is:
 *
 *		C = model_standard_mex(time, VIF, vp, kt, kep);
 *
 * Compilation:
 *
 *  mex -R2018a model_standard_mex.c
 *  or:
 *  1) uncomment the compiler directive MATLAB2015
 *  2) mex COPTIMFLAGS="\$COPTIMFLAGS -std=c99" model_standard_mex.c
 *
 * Yannick 2020
 *
 *========================================================*/

#include "mex.h"
#include <math.h>

#ifdef __GNU__
    #include <omp.h>
#endif

#ifndef MAXCORES
  #define MAXCORES 1
#endif 

// #define MATLAB2015

float dce_standard_value_float(
        float vp,   // vp
        float kt,  //Ktrans
        float kep,  // Kep
        int const point_index,   // time points to evaluate teh model at
        float const * T,            // time point vector
        float const * Cp            // VIF
        )
{
    // integral/convolution
    float convFunc = 0;
    for (int i = 1; i <= point_index; i++) {
        // simple sum rule
        float spacing = T[i] - T[i - 1];
        float Ct = Cp[i] * exp(- kep * (T[point_index]-T[i]));
        convFunc += (Ct * spacing);
    }
    
    float function_value = kt * convFunc + vp * Cp[point_index];
    return function_value;
}

double dce_standard_value_double(
        double vp,   // vp
        double kt,  //Ktrans
        double kep,  // Kep
        int const point_index,   // time points to evaluate teh model at
        double const * T,            // time point vector
        double const * Cp            // VIF
        )
{
    // integral/convolution
    double convFunc = 0;
    for (int i = 1; i <= point_index; i++) {
        // simple sum rule
        double spacing = T[i] - T[i - 1];
        double Ct = Cp[i] * exp(- kep * (T[point_index]-T[i]));
        convFunc += (Ct * spacing);
    }
    
    double function_value = kt * convFunc + vp * Cp[point_index];
    return function_value;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    size_t N, Nt, Na;
    mxClassID precision;
    float *Cpf, *timef, *vpf, *ktf, *kepf, *Cf;
    double *Cpd, *timed, *vpd, *ktd, *kepd, *Cd;
    
    /* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("STARDCE:model_standard:nrhs","Five inputs required: time, VIF, vp, kt, kep.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("STARDCE:model_standard:nlhs","One output required.");
    }
    
    // read all inputs
    Nt = mxGetN(prhs[0]);
    N = mxGetN(prhs[2]);
    if (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS) {
        precision = mxDOUBLE_CLASS;
        
#ifdef MATLAB2015
        // try to be backward compatible
        timed = mxGetPr(prhs[0]);
        Cpd = mxGetPr(prhs[1]);
        
        vpd = mxGetPr(prhs[2]);
        ktd = mxGetPr(prhs[3]);
        kepd = mxGetPr(prhs[4]);
#else
        timed = mxGetDoubles(prhs[0]);
        Cpd = mxGetDoubles(prhs[1]);
        
        vpd = mxGetDoubles(prhs[2]);
        ktd = mxGetDoubles(prhs[3]);
        kepd = mxGetDoubles(prhs[4]);
#endif
    }
    else {
        precision = mxSINGLE_CLASS;
#ifdef MATLAB2015
        // try to be backward compatible
        timef = mxGetData(prhs[0]);
        Cpf = mxGetData(prhs[1]);
        
        vpf = mxGetData(prhs[2]);
        ktf = mxGetData(prhs[3]);
        kepf = mxGetData(prhs[4]);
#else
        timef = mxGetSingles(prhs[0]);
        Cpf = mxGetSingles(prhs[1]);
        
        vpf = mxGetSingles(prhs[2]);
        ktf = mxGetSingles(prhs[3]);
        kepf = mxGetSingles(prhs[4]);
#endif
    }
    
    // output concentration
    plhs[0] = mxCreateNumericMatrix(N, Nt, precision, mxREAL);
    
    
#ifdef __GNU__
    /* Set number of threads */
    omp_set_num_threads(MAXCORES);
#endif
    
    
    if (precision == mxDOUBLE_CLASS){
#ifdef MATLAB2015
        Cd = mxGetPr(plhs[0]);
#else
        Cd = mxGetDoubles(plhs[0]);
#endif
        
        #pragma omp parallel for private(n,t)
        for (int n = 0; n < N; n++) {
            for (int t = 0; t < Nt; t++) {
                Cd[t * N + n] = dce_standard_value_double(vpd[n], ktd[n], kepd[n], t, timed, Cpd);
            }
        }
    } else {
#ifdef MATLAB2015
        Cf = mxGetData(plhs[0]);
#else
        Cf = mxGetSingles(plhs[0]);
#endif
        #pragma omp parallel for private(n,t)
        for (int n = 0; n < N; n++) {
            for (int t = 0; t < Nt; t++) {
                Cf[t * N + n] = dce_standard_value_float(vpf[n], ktf[n], kepf[n], t, timef, Cpf);
            }
        }
    }


    
}
