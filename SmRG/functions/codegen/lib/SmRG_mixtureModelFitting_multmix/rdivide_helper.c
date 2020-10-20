/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rdivide_helper.c
 *
 * Code generation for function 'rdivide_helper'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "rdivide_helper.h"

/* Function Definitions */
void rdivide_helper(const double x_data[], const int x_size[2], const double
                    y_data[], double z_data[], int z_size[2])
{
  int loop_ub;
  int i0;
  z_size[0] = 1;
  z_size[1] = x_size[1];
  loop_ub = x_size[0] * x_size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    z_data[i0] = x_data[i0] / y_data[i0];
  }
}

/* End of code generation (rdivide_helper.c) */
