/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_mixtureModelFitting_newPost.h
 *
 * Code generation for function 'SmRG_mixtureModelFitting_newPost'
 *
 */

#ifndef SMRG_MIXTUREMODELFITTING_NEWPOST_H
#define SMRG_MIXTUREMODELFITTING_NEWPOST_H

/* Include files */
#include "SmRG_mixtureModelFitting_newPost_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  extern void SmRG_mixtureModelFitting_newPost(const emxArray_real_T *Vin,
    double *K0_double, double *vB_double, double *mu_sk, double *rk, double tol,
    emxArray_real_T *p_tot, emxArray_real_T *a1);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (SmRG_mixtureModelFitting_newPost.h) */
