/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting.h"
#include "sum.h"
#include "combineVectorElements.h"

/* Function Definitions */
void b_sum(const emxArray_real_T *x, double y_data[], int y_size[2])
{
  int vlen;
  int k;
  vlen = x->size[0];
  if (x->size[0] == 0) {
    y_size[0] = 1;
    y_size[1] = 1;
    y_data[0] = 0.0;
  } else {
    y_size[0] = 1;
    y_size[1] = 1;
    y_data[0] = x->data[0];
    for (k = 2; k <= vlen; k++) {
      y_data[0] += x->data[k - 1];
    }
  }
}

double sum(const emxArray_real_T *x)
{
  return combineVectorElements(x);
}

/* End of code generation (sum.c) */
