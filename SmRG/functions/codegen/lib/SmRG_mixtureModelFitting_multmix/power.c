/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * power.c
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "power.h"
#include "SmRG_mixtureModelFitting_multmix_emxutil.h"

/* Function Definitions */
void power(const emxArray_real_T *a, emxArray_real_T *y)
{
  int nx;
  int k;
  nx = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = a->size[1];
  emxEnsureCapacity_real_T(y, nx);
  nx = a->size[1];
  for (k = 0; k < nx; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

/* End of code generation (power.c) */
