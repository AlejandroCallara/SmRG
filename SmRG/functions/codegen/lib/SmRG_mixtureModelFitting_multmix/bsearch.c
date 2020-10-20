/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * bsearch.c
 *
 * Code generation for function 'bsearch'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "bsearch.h"

/* Function Definitions */
int b_bsearch(const double x_data[], const int x_size[1], double xi)
{
  int n;
  int high_i;
  int low_ip1;
  int mid_i;
  high_i = x_size[0];
  n = 1;
  low_ip1 = 2;
  while (high_i > low_ip1) {
    mid_i = (n >> 1) + (high_i >> 1);
    if (((n & 1) == 1) && ((high_i & 1) == 1)) {
      mid_i++;
    }

    if (xi >= x_data[mid_i - 1]) {
      n = mid_i;
      low_ip1 = mid_i + 1;
    } else {
      high_i = mid_i;
    }
  }

  return n;
}

/* End of code generation (bsearch.c) */
