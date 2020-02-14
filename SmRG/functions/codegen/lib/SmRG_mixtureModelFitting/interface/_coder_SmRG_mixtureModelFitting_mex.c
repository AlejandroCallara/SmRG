/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_SmRG_mixtureModelFitting_mex.c
 *
 * Code generation for function '_coder_SmRG_mixtureModelFitting_mex'
 *
 */

/* Include files */
#include "_coder_SmRG_mixtureModelFitting_api.h"
#include "_coder_SmRG_mixtureModelFitting_mex.h"

/* Function Declarations */
static void c_SmRG_mixtureModelFitting_mexF(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[5]);

/* Function Definitions */
static void c_SmRG_mixtureModelFitting_mexF(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[5])
{
  const mxArray *outputs[6];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 5, 4,
                        24, "SmRG_mixtureModelFitting");
  }

  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 24,
                        "SmRG_mixtureModelFitting");
  }

  /* Call the function. */
  SmRG_mixtureModelFitting_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(SmRG_mixtureModelFitting_atexit);

  /* Module initialization. */
  SmRG_mixtureModelFitting_initialize();

  /* Dispatch the entry-point. */
  c_SmRG_mixtureModelFitting_mexF(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  SmRG_mixtureModelFitting_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_SmRG_mixtureModelFitting_mex.c) */
