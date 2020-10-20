/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "repmat.h"

/* Function Definitions */
void repmat(const cell_wrap_5 a[1], double varargin_1, cell_wrap_5 b_data[], int
            b_size[1])
{
  int i3;
  int itilerow;
  int loop_ub;
  int i4;
  i3 = (int)varargin_1;
  b_size[0] = (signed char)i3;
  for (itilerow = 0; itilerow < i3; itilerow++) {
    loop_ub = a[0].f1.size[0];
    b_data[itilerow].f1.size[0] = a[0].f1.size[0];
    for (i4 = 0; i4 < loop_ub; i4++) {
      b_data[itilerow].f1.data[i4] = a[0].f1.data[i4];
    }
  }
}

/* End of code generation (repmat.c) */
