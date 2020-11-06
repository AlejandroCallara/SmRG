/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_SmRG_mixtureModelFitting_newPost_api.h
 *
 * Code generation for function 'SmRG_mixtureModelFitting_newPost'
 *
 */

#ifndef _CODER_SMRG_MIXTUREMODELFITTING_NEWPOST_API_H
#define _CODER_SMRG_MIXTUREMODELFITTING_NEWPOST_API_H

/* Include files */
#include "emlrt.h"
#include "tmwtypes.h"
#include <string.h>

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  void SmRG_mixtureModelFitting_newPost(emxArray_real_T *Vin, real_T *K0_double,
    real_T *vB_double, real_T *mu_sk, real_T *rk, real_T tol, emxArray_real_T
    *p_tot, emxArray_real_T *a1);
  void SmRG_mixtureModelFitting_newPost_atexit(void);
  void SmRG_mixtureModelFitting_newPost_initialize(void);
  void SmRG_mixtureModelFitting_newPost_terminate(void);
  void SmRG_mixtureModelFitting_newPost_xil_shutdown(void);
  void SmRG_mixtureModelFitting_newPost_xil_terminate(void);
  void c_SmRG_mixtureModelFitting_newP(const mxArray * const prhs[6], int32_T
    nlhs, const mxArray *plhs[6]);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (_coder_SmRG_mixtureModelFitting_newPost_api.h) */
