/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mean.c
 *
 * Code generation for function 'mean'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "mean.h"
#include "combineVectorElements.h"

/* Function Definitions */
void b_mean(const emxArray_real_T *x, double y_data[], int y_size[2])
{
  int i1;
  int b_x;
  int loop_ub;
  c_combineVectorElements(x, y_data, y_size);
  i1 = y_size[0] * y_size[1];
  y_size[0] = 1;
  b_x = x->size[0];
  loop_ub = i1 - 1;
  for (i1 = 0; i1 <= loop_ub; i1++) {
    y_data[i1] /= (double)b_x;
  }
}

double c_mean(const emxArray_real_T *x)
{
  return b_combineVectorElements(x) / (double)x->size[0];
}

double d_mean(const double x_data[], const int x_size[2])
{
  double y;
  int vlen;
  int k;
  vlen = x_size[1];
  if (x_size[1] == 0) {
    y = 0.0;
  } else {
    y = x_data[0];
    for (k = 2; k <= vlen; k++) {
      y += x_data[k - 1];
    }
  }

  y /= (double)x_size[1];
  return y;
}

double mean(const emxArray_real_T *x)
{
  return combineVectorElements(x) / (double)x->size[1];
}

/* End of code generation (mean.c) */
