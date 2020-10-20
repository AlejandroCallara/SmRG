/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * log.c
 *
 * Code generation for function 'log'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "log.h"

/* Function Definitions */
void b_log(double *x)
{
  *x = log(*x);
}

void c_log(double x_data[], int x_size[2])
{
  int nx;
  int k;
  nx = x_size[1];
  for (k = 0; k < nx; k++) {
    x_data[k] = log(x_data[k]);
  }
}

/* End of code generation (log.c) */
