/*==========================================================
 * model_2cxm_mex.c - STARDCE toolbox
 *
 * Implements the DCE 2-compartment exchange model
 *
 *  Note: this is not really cpp but clang c-compiler is was not compatible with matlab after the recent MacOS update ... fix in the future! 
 *
 * The calling syntax is:
 *
 *		C = model_2cxm_mex(time, VIF, vp, ve, kt, fp );
 *
 * Compilation:
 *
 *  mex -R2018a model_2cxm_mex.c
 *  or:
 *  mex COPTIMFLAGS="\$COPTIMFLAGS -std=c99" model_2cxm_mex.c
 *
 * Yannick 2020
 * Copied from GPUfit by Sam Barnes
 *
 *========================================================*/

#include "mex.h"
#include <math.h>

float dce_2cxm_value_float (
        float vp,   // vp
        float ve, //ve
        float kt,  //Ktrans
        float fp,  // Fp
        int const point_index,   // time points to evaluate teh model at
        float const * T,            // time point vector
        float const * Cp            // VIF
        )
{
    // integral/convolution
    float PS;
    if(kt>=fp) {
        PS = 10e8;
    } else {
        PS = fp / ((fp / kt) - 1);
    }
    
    float convFunc = 0;
    float Tp = vp / (PS + fp);
    float Te = ve / PS;
    float Tb = vp / fp;
    float Kpos = 0.5 * (1/Tp + 1/Te + sqrt(pow(1/Tp + 1/Te,2) - 4 * 1/Te * 1/Tb));
    float Kneg = 0.5 * (1/Tp + 1/Te - sqrt(pow(1/Tp + 1/Te,2) - 4 * 1/Te * 1/Tb));
    float Eneg = (Kpos - 1/Tb) / (Kpos - Kneg);
    for (int i = 1; i < point_index; i++) {
        float spacing = T[i] - T[i - 1];
        float Ct =     Cp[i]     * (exp(-(T[point_index] - T[i])   * Kpos) + Eneg * (exp(-(T[point_index] - T[i])   * Kneg) - exp(-Kpos)));//(p2 * exp(-(T[point_index] - T[i])/Tp) + p0 * (1 - exp(-(T[point_index] - T[i])/Tp)));
        float Ctprev = Cp[i - 1] * (exp(-(T[point_index] - T[i-1]) * Kpos) + Eneg * (exp(-(T[point_index] - T[i-1]) * Kneg) - exp(-Kpos))); //(p2 * exp(-(T[point_index] - T[i-1])/Tp) + p0 * (1 - exp(-(T[point_index] - T[i-1])/Tp)));
        convFunc += ((Ct + Ctprev) / 2 * spacing);
    }
    float function_value = fp * convFunc;
    return function_value;
}

float dce_2cxm_value_double (
        double vp,   // vp
        double ve,  //ve
        double kt,  //Ktrans
        double fp,  // Fp
        int const point_index,   // time points to evaluate teh model at
        double const * T,            // time point vector
        double const * Cp            // VIF
        )
{
    // integral/convolution
    double PS;
    if(kt>=fp) {
        PS = 10e8;
    } else {
        PS = fp / ((fp / kt) - 1);
    }
    
    double convFunc = 0;
    double Tp = vp / (PS + fp);
    double Te = ve / PS;
    double Tb = vp / fp;
    double Kpos = 0.5 * (1/Tp + 1/Te + sqrt(pow(1/Tp + 1/Te,2) - 4 * 1/Te * 1/Tb));
    double Kneg = 0.5 * (1/Tp + 1/Te - sqrt(pow(1/Tp + 1/Te,2) - 4 * 1/Te * 1/Tb));
    double Eneg = (Kpos - 1/Tb) / (Kpos - Kneg);
    for (int i = 1; i < point_index; i++) {
        double spacing = T[i] - T[i - 1];
        double Ct =     Cp[i]     * (exp(-(T[point_index] - T[i])   * Kpos) + Eneg * (exp(-(T[point_index] - T[i])   * Kneg) - exp(-Kpos)));//(p2 * exp(-(T[point_index] - T[i])/Tp) + p0 * (1 - exp(-(T[point_index] - T[i])/Tp)));
        double Ctprev = Cp[i - 1] * (exp(-(T[point_index] - T[i-1]) * Kpos) + Eneg * (exp(-(T[point_index] - T[i-1]) * Kneg) - exp(-Kpos))); //(p2 * exp(-(T[point_index] - T[i-1])/Tp) + p0 * (1 - exp(-(T[point_index] - T[i-1])/Tp)));
        convFunc += ((Ct + Ctprev) / 2 * spacing);
    }
    double function_value = fp * convFunc;
    return function_value;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    size_t N, Nt;
    mxClassID precision;
    float *Cpf, *timef, *vpf, *vef, *ktf, *fpf, *Cf;
    double *Cpd, *timed, *vpd, *ved, *ktd, *fpd, *Cd;
    
    /* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("STARDCE:model_standard:nrhs","Six inputs required: time, VIF, vp, ve, kt, fp.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("STARDCE:model_standard:nlhs","One output required.");
    }
    
    // read all inputs
    Nt = mxGetN(prhs[0]);
    N = mxGetN(prhs[2]);
    if (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS) {
        precision = mxDOUBLE_CLASS;
        
        timed = mxGetDoubles(prhs[0]);
        Cpd = mxGetDoubles(prhs[1]);
        
        vpd = mxGetDoubles(prhs[2]);
        ved = mxGetDoubles(prhs[3]);
        ktd = mxGetDoubles(prhs[4]);
        fpd = mxGetDoubles(prhs[5]);
        
        // try to be backward compatible
//         timed = mxGetPr(prhs[0]);
//         Cpd = mxGetPr(prhs[1]);
//         
//         vpd = mxGetPr(prhs[2]);
//         ved = mxGetPr(prhs[3]);
//         ktd = mxGetPr(prhs[4]);
//         fpd = mxGetPr(prhs[5]);
    }
    else {
        precision = mxSINGLE_CLASS;
        
        timef = mxGetSingles(prhs[0]);
        Cpf = mxGetSingles(prhs[1]);
        
        vpf = mxGetSingles(prhs[2]);
        vef = mxGetSingles(prhs[3]);
        ktf = mxGetSingles(prhs[4]);
        fpf = mxGetSingles(prhs[5]);
        
        // try to be backward compatible
//         timef = mxGetData(prhs[0]);
//         Cpf = mxGetData(prhs[1]);
//         
//         vpf = mxGetData(prhs[2]);
//         vef = mxGetData(prhs[3]);
//         ktf = mxGetData(prhs[4]);
//         fpf = mxGetData(prhs[5]);
    }
    
    // output concentration
    plhs[0] = mxCreateNumericMatrix(N, Nt, precision, mxREAL);
    if (precision == mxDOUBLE_CLASS){
        Cd = mxGetDoubles(plhs[0]);
//         Cd = mxGetPr(plhs[0]);
        
        for (int n = 0; n < N; n++) {
            for (int t = 0; t < Nt; t++) {
                Cd[t * N + n] = dce_2cxm_value_double(vpd[n], vef[n], ktd[n], fpd[n], t, timed, Cpd);
            }
        }
    } else {
        Cf = mxGetSingles(plhs[0]);
//         Cf = mxGetData(plhs[0]);
        
        for (int n = 0; n < N; n++) {
            for (int t = 0; t < Nt; t++) {
                Cf[t * N + n] = dce_2cxm_value_float(vpf[n], vef[n], ktf[n], fpf[n], t, timef, Cpf);
            }
        }
    }


    
}
