/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_SmRG_mixtureModelFitting_newPost_mex.c
 *
 * Code generation for function 'SmRG_mixtureModelFitting_newPost'
 *
 */

/* Include files */
#include "_coder_SmRG_mixtureModelFitting_newPost_mex.h"
#include "_coder_SmRG_mixtureModelFitting_newPost_api.h"

/* Function Definitions */
void SmRG_mixtureModelFitting_newPost_mexFunction(int32_T nlhs, mxArray *plhs[6],
  int32_T nrhs, const mxArray *prhs[6])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *outputs[6];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 6, 4,
                        32, "SmRG_mixtureModelFitting_newPost");
  }

  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 32,
                        "SmRG_mixtureModelFitting_newPost");
  }

  /* Call the function. */
  c_SmRG_mixtureModelFitting_newP(prhs, nlhs, outputs);

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
  mexAtExit(&SmRG_mixtureModelFitting_newPost_atexit);

  /* Module initialization. */
  SmRG_mixtureModelFitting_newPost_initialize();

  /* Dispatch the entry-point. */
  SmRG_mixtureModelFitting_newPost_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  SmRG_mixtureModelFitting_newPost_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_SmRG_mixtureModelFitting_newPost_mex.c) */
