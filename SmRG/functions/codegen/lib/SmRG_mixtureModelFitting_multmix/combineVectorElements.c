/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * combineVectorElements.c
 *
 * Code generation for function 'combineVectorElements'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "combineVectorElements.h"

/* Function Definitions */
double b_combineVectorElements(const emxArray_real_T *x)
{
  double y;
  int vlen;
  int k;
  vlen = x->size[0];
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= vlen; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

void c_combineVectorElements(const emxArray_real_T *x, double y_data[], int
  y_size[2])
{
  int vlen;
  int npages;
  signed char sz_idx_1;
  int i;
  int xpageoffset;
  int k;
  vlen = x->size[0];
  if (x->size[0] == 0) {
    sz_idx_1 = (signed char)x->size[1];
    y_size[0] = 1;
    y_size[1] = sz_idx_1;
    if (0 <= sz_idx_1 - 1) {
      memset(&y_data[0], 0, (unsigned int)(sz_idx_1 * (int)sizeof(double)));
    }
  } else {
    npages = x->size[1];
    y_size[0] = 1;
    y_size[1] = x->size[1];
    for (i = 0; i < npages; i++) {
      xpageoffset = i * x->size[0];
      y_data[i] = x->data[xpageoffset];
      for (k = 2; k <= vlen; k++) {
        y_data[i] += x->data[(xpageoffset + k) - 1];
      }
    }
  }
}

double combineVectorElements(const emxArray_real_T *x)
{
  double y;
  int vlen;
  int k;
  vlen = x->size[1];
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= vlen; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

/* End of code generation (combineVectorElements.c) */
