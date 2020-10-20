/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nbinlike_mu.c
 *
 * Code generation for function 'nbinlike_mu'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "nbinlike_mu.h"
#include "gammaln.h"

/* Function Declarations */
static double rt_roundd_snf(double u);

/* Function Definitions */
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

double nbinlike_mu(const double params[2], double data)
{
  double logL;
  double r;
  double p;
  double d0;
  double d1;
  double d2;

  /* |==================================================================================== */
  /* |NBINLIKE_MU Negative of the negative binomial log-likelihood function. */
  /* | */
  /* |   L = NBINLIKE_MU(PARAMS,DATA) returns the negative of the negative binomial */
  /* |                                log-likelihood function for the parameters */
  /* |                                PARAMS(1) = \mu and PARAMS(2) = \alpha, given DATA. */
  /* | */
  /* |   [LOGL, AVAR] = NBINLIKE_MU(PARAMS,DATA) adds the inverse of Fisher's information  */
  /* |                                           matrix, AVAR. If the input parameter  */
  /* |                                           values in PARAMS are the maximum likelihood */
  /* |                                           estimates, the diagonal elements of AVAR */
  /* |                                           are their asymptotic variances.   */
  /* |                                           AVAR is based on the observed Fisher's */
  /* |                                           information, not the expected information. */
  /* | */
  /* |   NBINLIKE_MU is a utility function for maximum likelihood estimation.  */
  /* | */
  /* |  Last revision: */
  /* |  22 May 2018 */
  /* |  Michele Scipioni, University of Pisa */
  /* | */
  /* |==================================================================================== */
  /*  if nargin < 2,  */
  /*      error(message('stats:nbinlike:TooFewInputs'));  */
  /*  end */
  /*  if min(n,m) > 1 */
  /*      error(message('stats:nbinlike:InvalidData')); */
  /*  end */
  r = 1.0 / params[1];
  p = r / (r + params[0]);

  /*  Out of range or missing parameters or data return NaN.  Infinite */
  /*  values for R correspond to a Poisson, but its mean cannot be determined */
  /*  from the (R,P) parametrization. */
  if ((!(0.0 < r)) || (rtIsInf(r) || rtIsNaN(r)) || (!(0.0 < p)) || (!(p <= 1.0))
      || ((!(0.0 <= data)) || (rtIsInf(data) || rtIsNaN(data)) || (!(data ==
         rt_roundd_snf(data))))) {
    logL = rtNaN;
  } else {
    d0 = r + data;
    gammaln(&d0);
    d1 = data + 1.0;
    gammaln(&d1);
    d2 = r;
    gammaln(&d2);
    logL = -(((d0 - d1) - d2) + r * log(p)) - data * log(1.0 - p);
  }

  return logL;
}

/* End of code generation (nbinlike_mu.c) */
