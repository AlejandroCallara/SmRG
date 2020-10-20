/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_SmRG_mixtureModelFitting_newPost_api.h
 *
 * Code generation for function '_coder_SmRG_mixtureModelFitting_newPost_api'
 *
 */

#ifndef _CODER_SMRG_MIXTUREMODELFITTING_NEWPOST_API_H
#define _CODER_SMRG_MIXTUREMODELFITTING_NEWPOST_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_SmRG_mixtureModelFitting_newPost_api.h"

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

/* Function Declarations */
extern void SmRG_mixtureModelFitting_newPost(emxArray_real_T *Vin, real_T
  *K0_double, real_T *vB_double, real_T *mu_sk, real_T *rk, emxArray_real_T
  *p_tot, emxArray_real_T *a1);
extern void SmRG_mixtureModelFitting_newPost_api(const mxArray * const prhs[5],
  int32_T nlhs, const mxArray *plhs[6]);
extern void SmRG_mixtureModelFitting_newPost_atexit(void);
extern void SmRG_mixtureModelFitting_newPost_initialize(void);
extern void SmRG_mixtureModelFitting_newPost_terminate(void);
extern void SmRG_mixtureModelFitting_newPost_xil_shutdown(void);
extern void SmRG_mixtureModelFitting_newPost_xil_terminate(void);

#endif

/* End of code generation (_coder_SmRG_mixtureModelFitting_newPost_api.h) */
