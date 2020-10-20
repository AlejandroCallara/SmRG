/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_SmRG_mixtureModelFitting_multmix_api.c
 *
 * Code generation for function '_coder_SmRG_mixtureModelFitting_multmix_api'
 *
 */

/* Include files */
#include <string.h>
#include "tmwtypes.h"
#include "_coder_SmRG_mixtureModelFitting_multmix_api.h"
#include "_coder_SmRG_mixtureModelFitting_multmix_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131482U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "SmRG_mixtureModelFitting_multmix",  /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static const mxArray *b_emlrt_marshallOut(const real_T u);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *K0_double,
  const char_T *identifier);
static void c_emlrt_marshallOut(const real_T u_data[], const int32_T u_size[2],
  const mxArray *y);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *mu_sk, const
  char_T *identifier, real_T **y_data, int32_T y_size[2]);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *Vin, const
  char_T *identifier, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int32_T oldNumel);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T **y_data, int32_T y_size[2]);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T **ret_data, int32_T ret_size[2]);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  g_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m1;
  y = NULL;
  m1 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m1);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *K0_double,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(K0_double), &thisId);
  emlrtDestroyArray(&K0_double);
  return y;
}

static void c_emlrt_marshallOut(const real_T u_data[], const int32_T u_size[2],
  const mxArray *y)
{
  emlrtMxSetData((mxArray *)y, (void *)&u_data[0]);
  emlrtSetDimensions((mxArray *)y, u_size, 2);
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *mu_sk, const
  char_T *identifier, real_T **y_data, int32_T y_size[2])
{
  emlrtMsgIdentifier thisId;
  real_T *r0;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(mu_sk), &thisId, &r0, y_size);
  *y_data = r0;
  emlrtDestroyArray(&mu_sk);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *Vin, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(Vin), &thisId, y);
  emlrtDestroyArray(&Vin);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[2] = { 0, 0 };

  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m0, &u->data[0]);
  emlrtSetDimensions((mxArray *)m0, u->size, 2);
  emlrtAssign(&y, m0);
  return y;
}

static void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int32_T oldNumel)
{
  int32_T newNumel;
  int32_T i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = emlrtCallocMex((uint32_T)i, sizeof(real_T));
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, sizeof(real_T) * oldNumel);
      if (emxArray->canFreeData) {
        emlrtFreeMex(emxArray->data);
      }
    }

    emxArray->data = (real_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (real_T *)NULL) && (*pEmxArray)->canFreeData) {
      emlrtFreeMex((*pEmxArray)->data);
    }

    emlrtFreeMex((*pEmxArray)->size);
    emlrtFreeMex(*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)emlrtMallocMex(sizeof(emxArray_real_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2012b(sp, (void *)pEmxArray, (void (*)(void *))
      emxFree_real_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex(sizeof(int32_T) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T **y_data, int32_T y_size[2])
{
  real_T *r1;
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, &r1, y_size);
  *y_data = r1;
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { 1, -1 };

  const boolean_T bv0[2] = { false, true };

  int32_T iv1[2];
  int32_T i0;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv0[0],
    iv1);
  ret->allocatedSize = iv1[0] * iv1[1];
  i0 = ret->size[0] * ret->size[1];
  ret->size[0] = iv1[0];
  ret->size[1] = iv1[1];
  emxEnsureCapacity_real_T(ret, i0);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T **ret_data, int32_T ret_size[2])
{
  static const int32_T dims[2] = { 1, 4 };

  const boolean_T bv1[2] = { false, true };

  int32_T iv2[2];
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv1[0],
    iv2);
  ret_size[0] = iv2[0];
  ret_size[1] = iv2[1];
  *ret_data = (real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
}

void SmRG_mixtureModelFitting_multmix_api(const mxArray *prhs[5], int32_T nlhs,
  const mxArray *plhs[6])
{
  emxArray_real_T *Vin;
  emxArray_real_T *p_tot;
  emxArray_real_T *a1;
  const mxArray *prhs_copy_idx_0;
  const mxArray *prhs_copy_idx_1;
  const mxArray *prhs_copy_idx_2;
  const mxArray *prhs_copy_idx_3;
  const mxArray *prhs_copy_idx_4;
  real_T K0_double;
  real_T vB_double;
  real_T (*mu_sk_data)[4];
  int32_T mu_sk_size[2];
  real_T (*rk_data)[4];
  int32_T rk_size[2];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &Vin, 2, true);
  emxInit_real_T(&st, &p_tot, 2, true);
  emxInit_real_T(&st, &a1, 2, true);
  prhs_copy_idx_0 = prhs[0];
  prhs_copy_idx_1 = prhs[1];
  prhs_copy_idx_2 = prhs[2];
  prhs_copy_idx_3 = emlrtProtectR2012b(prhs[3], 3, true, 4);
  prhs_copy_idx_4 = emlrtProtectR2012b(prhs[4], 4, true, 4);

  /* Marshall function inputs */
  Vin->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_0), "Vin", Vin);
  K0_double = c_emlrt_marshallIn(&st, emlrtAliasP(prhs_copy_idx_1), "K0_double");
  vB_double = c_emlrt_marshallIn(&st, emlrtAliasP(prhs_copy_idx_2), "vB_double");
  e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_3), "mu_sk", (real_T **)
                     &mu_sk_data, mu_sk_size);
  e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_4), "rk", (real_T **)&rk_data,
                     rk_size);

  /* Invoke the target function */
  SmRG_mixtureModelFitting_multmix(Vin, &K0_double, &vB_double, *mu_sk_data,
    mu_sk_size, *rk_data, rk_size, p_tot, a1);

  /* Marshall function outputs */
  p_tot->canFreeData = false;
  plhs[0] = emlrt_marshallOut(p_tot);
  emxFree_real_T(&p_tot);
  emxFree_real_T(&Vin);
  if (nlhs > 1) {
    a1->canFreeData = false;
    plhs[1] = emlrt_marshallOut(a1);
  }

  emxFree_real_T(&a1);
  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(K0_double);
  }

  if (nlhs > 3) {
    plhs[3] = b_emlrt_marshallOut(vB_double);
  }

  if (nlhs > 4) {
    c_emlrt_marshallOut(*mu_sk_data, mu_sk_size, prhs_copy_idx_3);
    plhs[4] = prhs_copy_idx_3;
  }

  if (nlhs > 5) {
    c_emlrt_marshallOut(*rk_data, rk_size, prhs_copy_idx_4);
    plhs[5] = prhs_copy_idx_4;
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

void SmRG_mixtureModelFitting_multmix_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  SmRG_mixtureModelFitting_multmix_xil_terminate();
  SmRG_mixtureModelFitting_multmix_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void SmRG_mixtureModelFitting_multmix_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void SmRG_mixtureModelFitting_multmix_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_SmRG_mixtureModelFitting_multmix_api.c) */
