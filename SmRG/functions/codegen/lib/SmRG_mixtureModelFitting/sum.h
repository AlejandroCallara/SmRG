/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.h
 *
 * Code generation for function 'sum'
 *
 */

#ifndef SUM_H
#define SUM_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "SmRG_mixtureModelFitting_types.h"

/* Function Declarations */
extern void b_sum(const emxArray_real_T *x, double y_data[], int y_size[2]);
extern double sum(const emxArray_real_T *x);

#endif

/* End of code generation (sum.h) */
