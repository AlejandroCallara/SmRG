/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mean.h
 *
 * Code generation for function 'mean'
 *
 */

#ifndef MEAN_H
#define MEAN_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "SmRG_mixtureModelFitting_multmix_types.h"

/* Function Declarations */
extern void b_mean(const emxArray_real_T *x, double y_data[], int y_size[2]);
extern double c_mean(const emxArray_real_T *x);
extern double d_mean(const double x_data[], const int x_size[2]);
extern double mean(const emxArray_real_T *x);

#endif

/* End of code generation (mean.h) */
