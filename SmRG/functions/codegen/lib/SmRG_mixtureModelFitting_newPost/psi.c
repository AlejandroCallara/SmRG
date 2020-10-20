/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * psi.c
 *
 * Code generation for function 'psi'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_newPost.h"
#include "psi.h"

/* Function Definitions */
double minusDoubleScalarPsi(double x)
{
  double rval;
  double xdmln;
  double xdmy;
  int xinc;
  double tt;
  double tst;
  double ta;
  double s;
  double tk;
  int floorxinc;
  boolean_T exitg1;
  int i;
  double trm;
  static const double dv0[22] = { 1.0, -0.5, 0.16666666666666666,
    -0.033333333333333333, 0.023809523809523808, -0.033333333333333333,
    0.07575757575757576, -0.2531135531135531, 1.1666666666666667,
    -7.0921568627450977, 54.971177944862156, -529.12424242424242,
    6192.123188405797, -86580.253113553117, 1.4255171666666667E+6,
    -2.7298231067816094E+7, 6.015808739006424E+8, -1.5116315767092157E+10,
    4.2961464306116669E+11, -1.3711655205088332E+13, 4.8833231897359319E+14,
    -1.9296579341940068E+16 };

  if (rtIsNaN(x)) {
    rval = rtNaN;
  } else if (x == 0.0) {
    rval = rtInf;
  } else if (rtIsInf(x)) {
    rval = rtMinusInf;
  } else if (x < 2.2204460492503131E-16) {
    rval = 1.0 / x;
  } else {
    xdmln = log(x);
    xdmy = x;
    xinc = 0;
    if (x < 9.0) {
      xdmy = floor(x);
      xinc = 9 - (int)xdmy;
      xdmy = x + (9.0 - xdmy);
      xdmln = log(xdmy);
    }

    tt = 0.5 / xdmy;
    tst = 2.2204460492503131E-16 * tt;
    xdmy = 1.0 / (xdmy * xdmy);
    ta = 0.5 * xdmy;
    s = ta / 6.0;
    if (s >= tst) {
      tk = 2.0;
      floorxinc = 0;
      exitg1 = false;
      while ((!exitg1) && (floorxinc < 19)) {
        ta = ta * ((tk + 1.0) / (tk + 1.0)) * (tk / (tk + 2.0)) * xdmy;
        trm = ta * dv0[floorxinc + 3];
        if (fabs(trm) < tst) {
          exitg1 = true;
        } else {
          s += trm;
          tk += 2.0;
          floorxinc++;
        }
      }
    }

    s = (s + tt) * exp(-(0.0 * xdmln));
    if (xinc > 0) {
      floorxinc = xinc - 1;
      for (i = 0; i <= floorxinc; i++) {
        s += 1.0 / (x + (double)((xinc - i) - 1));
      }

      rval = s - xdmln;
    } else {
      rval = s - xdmln;
    }
  }

  return -rval;
}

/* End of code generation (psi.c) */
