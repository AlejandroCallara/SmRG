/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_mixtureModelFitting.h
 *
 * Code generation for function 'SmRG_mixtureModelFitting'
 *
 */

#ifndef SMRG_MIXTUREMODELFITTING_H
#define SMRG_MIXTUREMODELFITTING_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "SmRG_mixtureModelFitting_types.h"

/* Function Declarations */
extern void SmRG_mixtureModelFitting(const emxArray_real_T *Vin, double
  *K0_double, const double *vB_double, double *mu_sk, double *rk,
  emxArray_real_T *p_tot, emxArray_real_T *a1);

#endif

/* End of code generation (SmRG_mixtureModelFitting.h) */
