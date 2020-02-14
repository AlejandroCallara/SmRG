/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_mixtureModelFitting.c
 *
 * Code generation for function 'SmRG_mixtureModelFitting'
 *
 */

/* Include files */
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting.h"
#include "SmRG_mixtureModelFitting_emxutil.h"
#include "combineVectorElements.h"
#include "mean.h"
#include "sum.h"
#include "power.h"
#include "exp.h"
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

void SmRG_mixtureModelFitting(const emxArray_real_T *Vin, double *K0_double,
  const double *vB_double, double *mu_sk, double *rk, emxArray_real_T *p_tot,
  emxArray_real_T *a1)
{
  int i0;
  int loop_ub;
  emxArray_uint32_T *i_ok;
  unsigned int u0;
  double pb;
  double pa;
  int varargin_1;
  int n_of_pixel;
  double log_lik_a[5000];
  double log_lik_b[5000];
  int count;
  emxArray_real_T *psi_a;
  emxArray_real_T *psi_b;
  emxArray_real_T *b;
  emxArray_real_T *a;
  emxArray_real_T *b_b;
  emxArray_real_T *b_Vin;
  emxArray_real_T *b_a1;
  boolean_T exitg1;
  emxArray_int32_T *r0;
  double r;
  double B_tmp_data[1];
  int B_tmp_size[2];
  double z;
  double nlogL_tmp;
  double d0;
  double d1;
  double d2;

  /*  SmRG_mixtureModelFitting: */
  /*            fits model described in [1] on data input Vin. The model */
  /*            describes a single class k of signal pixels (negative-binomial) */
  /*            and the background (gaussian).For more details refer to [1]. */
  /*  */
  /*  Syntax: */
  /*            [p_tot,a1]=SmRG_mixtureModelFitting(vin) */
  /*            [p_tot,a1,K0_double,vB,mu_sk,rk]=SmRG_mixtureModelFitting(vin) */
  /*  Inputs: */
  /*            Vin: 2D image (type: double)* */
  /*  */
  /*  Outputs: */
  /*            p_tot:  posterior probabilities for each element in Vin(:) of */
  /*                    belonging to neg-bin distribution. The posterior */
  /*                    probabilities of belonging to the gaussian distribution */
  /*                    are simply: */
  /*                                p_gaussian = 1-p_tot */
  /*            a1:     posterior probabilities of voxels belonging to the */
  /*                    negbin distribution */
  /*            K0_double: Gaussian mean */
  /*            vB:     Gaussian variance */
  /*            mu_sk:  negbin mean */
  /*            rk:     negbin dispersion parameter */
  /*  */
  /*  */
  /*            * note: the function allows any N-dimensional input. When */
  /*            called from SmRG_regionGrowing Vin is a 2D image. However, one */
  /*            can use this function to fit the mixture model of [1] on any */
  /*            kind of input, regardless of dimensionality. */
  /*  */
  /* [1]: Calapez,A. and Rosa,A. (2010) A statistical pixel intensity model */
  /*      for segmentation of confocal laser scanning microscopy images. */
  /*      IEEE Trans. Image Process., 19, 2408–2418. */
  /*  % check inputs */
  /*  if nargin<1 */
  /*      help SmRG_mixtureModelFitting */
  /*      return */
  /*  end */
  /*   */
  /*  % if nargin<2 */
  /*  %     % if background not specified never do check on background */
  /*  %     background = 0; */
  /*  % end */
  /*   */
  /*  if sum(Vin(:)<0) */
  /*      error ('only images with positive values are allowed') */
  /*  end */
  /*  check for nans in data */
  /*  [Vd,i_nan,i_ok]=SmRG_workWithNans(vin); */
  /*  pre-allocate p_tot */
  i0 = p_tot->size[0] * p_tot->size[1];
  p_tot->size[0] = 1;
  p_tot->size[1] = Vin->size[0] * Vin->size[1];
  emxEnsureCapacity_real_T(p_tot, i0);
  loop_ub = Vin->size[0] * Vin->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    p_tot->data[i0] = 0.0;
  }

  i0 = a1->size[0] * a1->size[1];
  a1->size[0] = 1;
  a1->size[1] = Vin->size[0] * Vin->size[1];
  emxEnsureCapacity_real_T(a1, i0);
  loop_ub = Vin->size[0] * Vin->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    a1->data[i0] = 0.0;
  }

  emxInit_uint32_T(&i_ok, 2);
  if (Vin->size[0] * Vin->size[1] < 1) {
    i_ok->size[0] = 1;
    i_ok->size[1] = 0;
  } else {
    u0 = (unsigned int)(Vin->size[0] * Vin->size[1]);
    i0 = i_ok->size[0] * i_ok->size[1];
    i_ok->size[0] = 1;
    i_ok->size[1] = (int)u0;
    emxEnsureCapacity_uint32_T(i_ok, i0);
    loop_ub = (int)u0 - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      i_ok->data[i0] = 1U + i0;
    }
  }

  /*  if exist('vB') */
  /*      vB_double = double(vB); */
  /*  end */
  /*  Mixing prior probability (0.5 0.5) */
  pb = 0.5;
  pa = 0.5;

  /*  prior */
  /*  EM initial conditions */
  varargin_1 = Vin->size[0] * Vin->size[1] - 1;
  n_of_pixel = Vin->size[0] * Vin->size[1];
  memset(&log_lik_a[0], 0, 5000U * sizeof(double));
  memset(&log_lik_b[0], 0, 5000U * sizeof(double));
  log_lik_a[0] = 0.0;
  log_lik_b[0] = 0.0;
  log_lik_a[1] = 10.0;
  log_lik_b[1] = 10.0;
  count = 1;

  /*  figure */
  /*      h = histogram(Vin) */
  /*      hold on; */
  /*      asse_x =1:1000;% linspace(h.BinLimits(1),h.BinLimits(2),2000); */
  /*      hold on */
  /*      plot(asse_x,90000*normpdf(asse_x,K0_double,sqrt(vB))); */
  /*      plot(asse_x,90000*nbinpdf_mu(asse_x,mu_sk,1/rk)); */
  /*      drawnow */
  /*  mixture, EM */
  /*  iterates until convergence */
  emxInit_real_T(&psi_a, 1);
  emxInit_real_T(&psi_b, 1);
  emxInit_real_T(&b, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&b_b, 1);
  emxInit_real_T(&b_Vin, 1);
  emxInit_real_T(&b_a1, 2);
  exitg1 = false;
  while ((!exitg1) && ((fabs(log_lik_a[count] - log_lik_a[count - 1]) > 1.0E-6) &&
                       (fabs(log_lik_b[count] - log_lik_b[count - 1]) > 1.0E-6)))
  {
    /*  initialize log likelihood */
    i0 = psi_a->size[0];
    psi_a->size[0] = varargin_1 + 1;
    emxEnsureCapacity_real_T(psi_a, i0);
    for (i0 = 0; i0 <= varargin_1; i0++) {
      psi_a->data[i0] = 0.0;
    }

    i0 = psi_b->size[0];
    psi_b->size[0] = n_of_pixel;
    emxEnsureCapacity_real_T(psi_b, i0);
    for (i0 = 0; i0 < n_of_pixel; i0++) {
      psi_b->data[i0] = 0.0;
    }

    /*  initialize log posterior */
    i0 = b->size[0] * b->size[1];
    b->size[0] = n_of_pixel;
    b->size[1] = 1;
    emxEnsureCapacity_real_T(b, i0);
    for (i0 = 0; i0 < n_of_pixel; i0++) {
      b->data[i0] = 0.0;
    }

    i0 = a->size[0] * a->size[1];
    a->size[0] = n_of_pixel;
    a->size[1] = 1;
    emxEnsureCapacity_real_T(a, i0);
    for (i0 = 0; i0 < n_of_pixel; i0++) {
      a->data[i0] = 0.0;
    }

    /*  E-step */
    for (loop_ub = 0; loop_ub <= varargin_1; loop_ub++) {
      /*  cast */
      /*          K0_double = double(K0); vB_double = double(vB); */
      /*   Gaussian log likelihood */
      r = sqrt(*vB_double);

      /* NORMLIKE Negative log-likelihood for the normal distribution. */
      /*    NLOGL = NORMLIKE(PARAMS,DATA) returns the negative of the log-likelihood */
      /*    for the normal distribution, evaluated at parameters PARAMS(1) = MU and */
      /*    PARAMS(2) = SIGMA, given DATA.  NLOGL is a scalar. */
      /*  */
      /*    [NLOGL, AVAR] = NORMLIKE(PARAMS,DATA) returns the inverse of Fisher's */
      /*    information matrix, AVAR.  If the input parameter values in PARAMS are */
      /*    the maximum likelihood estimates, the diagonal elements of AVAR are */
      /*    their asymptotic variances.  AVAR is based on the observed Fisher's */
      /*    information, not the expected information. */
      /*  */
      /*    [...] = NORMLIKE(PARAMS,DATA,CENSORING) accepts a boolean vector of the */
      /*    same size as DATA that is 1 for observations that are right-censored */
      /*    and 0 for observations that are observed exactly. */
      /*  */
      /*    [...] = NORMLIKE(PARAMS,DATA,CENSORING,FREQ) accepts a frequency vector */
      /*    of the same size as DATA.  FREQ typically contains integer frequencies */
      /*    for the corresponding elements in DATA, but may contain any non-integer */
      /*    non-negative values.  Pass in [] for CENSORING to use its default */
      /*    value. */
      /*  */
      /*    See also NORMCDF, NORMFIT, NORMINV, NORMPDF, NORMRND, NORMSTAT. */
      /*    References: */
      /*       [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical */
      /*           Distributions, 2nd ed., Wiley, 170pp. */
      /*       [2] Lawless, J.F. (1982) Statistical Models and Methods for Lifetime */
      /*           Data, Wiley, New York, 580pp. */
      /*       [3} Meeker, W.Q. and L.A. Escobar (1998) Statistical Methods for */
      /*           Reliability Data, Wiley, New York, 680pp. */
      /*    Copyright 1993-2014 The MathWorks, Inc. */
      /*  if nargin < 2 */
      /*      error(message('stats:normlike:TooFewInputs')); */
      /*  elseif numel(data) > length(data) */
      /*      error(message('stats:normlike:InvalidData')); */
      /*  end */
      /*  if numel(params)~=2 */
      /*      error(message('stats:probdists:WrongParameterLength',2)); */
      /*  end */
      /*  if nargin < 3 || isempty(censoring) */
      /*      censoring = 0; % make this a scalar, will expand when needed */
      /*  elseif ~isequal(size(data), size(censoring)) */
      /*      error(message('stats:normlike:InputSizeMismatchCensoring')); */
      /*  end */
      /*  if nargin < 4 || isempty(freq) */
      /*      freq = 1; % make this a scalar, will expand when needed */
      /*  elseif isequal(size(data), size(freq)) */
      /*      zerowgts = find(freq == 0); */
      /*      if numel(zerowgts) > 0 */
      /*          data(zerowgts) = []; */
      /*          if numel(censoring)==numel(freq), censoring(zerowgts) = []; end */
      /*          freq(zerowgts) = []; */
      /*      end */
      /*  else */
      /*      error(message('stats:normlike:InputSizeMismatchFreq')); */
      /*  end */
      /*  Return NaN for out of range parameters. */
      if (r <= 0.0) {
        r = rtNaN;
      }

      /*  Compute the individual log-likelihood terms. */
      z = (Vin->data[loop_ub] - *K0_double) / r;

      /*  Sum up the individual contributions, and return the negative log-likelihood. */
      nlogL_tmp = -0.5 * z * z - log(2.5066282746310002 * r);

      /*  Compute the negative hessian at the parameter values, and invert to get */
      /*  the observed information matrix. */
      /*  if nargout == 2 */
      /*      dL11 = -ones(size(z), 'like', z); */
      /*      dL12 = -2.*z; */
      /*      dL22 = 1 - 3.*z.*z; */
      /*      if ncen > 0 */
      /*          dlogScen = exp(-.5*zcen.*zcen) ./ (sqrt(2*pi) .* Scen); */
      /*          d2logScen = dlogScen .* (dlogScen - zcen); */
      /*          dL11(cens) = -d2logScen; */
      /*          dL12(cens) = -dlogScen - zcen.*d2logScen; */
      /*          dL22(cens) = -zcen .* (2.*dlogScen + zcen.*d2logScen); */
      /*      end */
      /*      nH11 = -sum(freq .* dL11); */
      /*      nH12 = -sum(freq .* dL12); */
      /*      nH22 = -sum(freq .* dL22); */
      /*      avar =  (sigma.^2) * [nH22 -nH12; -nH12 nH11] / (nH11*nH22 - nH12*nH12); */
      /*  end */
      /*          if isnan(psi_b_tmp);keyboard;end */
      /* if (Vtmp-K0_double)<0;K0_double=double(min(V));end */
      psi_b->data[loop_ub] = nlogL_tmp;

      /*   Neg-Bin log likelihood */
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
      r = 1.0 / (1.0 / *rk);
      z = r / (r + *mu_sk);

      /*  Out of range or missing parameters or data return NaN.  Infinite */
      /*  values for R correspond to a Poisson, but its mean cannot be determined */
      /*  from the (R,P) parametrization. */
      if ((!(0.0 < r)) || (rtIsInf(r) || rtIsNaN(r)) || (!(0.0 < z)) || (!(z <=
            1.0)) || ((!(0.0 <= Vin->data[loop_ub])) || (rtIsInf(Vin->
             data[loop_ub]) || rtIsNaN(Vin->data[loop_ub])) || (!(Vin->
             data[loop_ub] == rt_roundd_snf(Vin->data[loop_ub]))))) {
        r = rtNaN;
      } else {
        d0 = r + Vin->data[loop_ub];
        gammaln(&d0);
        d1 = Vin->data[loop_ub] + 1.0;
        gammaln(&d1);
        d2 = r;
        gammaln(&d2);
        r = -(((d0 - d1) - d2) + r * log(z)) - Vin->data[loop_ub] * log(1.0 - z);
      }

      /*          if isnan(psi_a_tmp);keyboard;end */
      psi_a->data[loop_ub] = -r;

      /*  log posterior */
      d0 = log(exp(nlogL_tmp) * pb + exp(-r) * pa);
      b->data[loop_ub] = (nlogL_tmp + log(pb)) - d0;
      a->data[loop_ub] = (-r + log(pa)) - d0;
    }

    b_exp(b);
    b_exp(a);
    i0 = a1->size[0] * a1->size[1];
    a1->size[0] = a->size[0];
    a1->size[1] = 1;
    emxEnsureCapacity_real_T(a1, i0);
    loop_ub = a->size[0] * a->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      a1->data[i0] = a->data[i0];
    }

    /*      if sum(isnan(b1))~=0;keyboard;end */
    /*      if sum(isnan(a1))~=0;keyboard;end */
    /*  M-step */
    /*  update Gaussian Parameters */
    b_sum(b, B_tmp_data, B_tmp_size);
    loop_ub = b->size[0];
    i0 = b_b->size[0];
    b_b->size[0] = loop_ub;
    emxEnsureCapacity_real_T(b_b, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_b->data[i0] = b->data[i0] * Vin->data[i0];
    }

    *K0_double = sum(b_b) / B_tmp_data[0];

    /*  update Neg-Bin Parameters */
    b_sum(a1, B_tmp_data, B_tmp_size);
    loop_ub = a1->size[0];
    i0 = b_b->size[0];
    b_b->size[0] = loop_ub;
    emxEnsureCapacity_real_T(b_b, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_b->data[i0] = a1->data[i0] * Vin->data[i0];
    }

    *mu_sk = sum(b_b) / B_tmp_data[0];
    loop_ub = Vin->size[0] * Vin->size[1];
    i0 = b_Vin->size[0];
    b_Vin->size[0] = loop_ub;
    emxEnsureCapacity_real_T(b_Vin, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_Vin->data[i0] = Vin->data[i0] - *mu_sk;
    }

    power(b_Vin, b_b);
    loop_ub = a1->size[0];
    i0 = b_b->size[0];
    b_b->size[0] = loop_ub;
    emxEnsureCapacity_real_T(b_b, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_b->data[i0] *= a1->data[i0];
    }

    r = sum(b_b) / B_tmp_data[0];

    /*    pk = mu_sk/v_sk; */
    /*  for numerical instability, to fix in future */
    if (r < *mu_sk) {
      r = 2.0 * *mu_sk;
    }

    *rk = *mu_sk * *mu_sk / (r - *mu_sk);

    /*  today's posterior becomes tomorrow's prior :) */
    loop_ub = a1->size[0];
    i0 = b_a1->size[0] * b_a1->size[1];
    b_a1->size[0] = 1;
    b_a1->size[1] = loop_ub;
    emxEnsureCapacity_real_T(b_a1, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_a1->data[i0] = a1->data[i0];
    }

    pa = mean(b_a1);
    loop_ub = b->size[0];
    i0 = b_a1->size[0] * b_a1->size[1];
    b_a1->size[0] = 1;
    b_a1->size[1] = loop_ub;
    emxEnsureCapacity_real_T(b_a1, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_a1->data[i0] = b->data[i0];
    }

    pb = mean(b_a1);

    /*  update for tollerance check */
    count++;
    log_lik_a[count] = combineVectorElements(psi_a) / (double)psi_a->size[0];
    log_lik_b[count] = combineVectorElements(psi_b) / (double)psi_b->size[0];

    /*  prevent infinite loop */
    if (count + 1 > 4999) {
      /*          keyboard */
      exitg1 = true;
    }
  }

  emxFree_real_T(&b_a1);
  emxFree_real_T(&b_Vin);
  emxFree_real_T(&b_b);
  emxFree_real_T(&a);
  emxFree_real_T(&b);
  emxFree_real_T(&psi_b);
  emxFree_real_T(&psi_a);
  emxInit_int32_T(&r0, 2);
  i0 = r0->size[0] * r0->size[1];
  r0->size[0] = 1;
  r0->size[1] = i_ok->size[1];
  emxEnsureCapacity_int32_T(r0, i0);
  loop_ub = i_ok->size[0] * i_ok->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    r0->data[i0] = (int)i_ok->data[i0];
  }

  emxFree_uint32_T(&i_ok);
  loop_ub = r0->size[0] * r0->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    p_tot->data[r0->data[i0] - 1] = a1->data[i0];
  }

  emxFree_int32_T(&r0);

  /*  p_tot(i_nan)=1; */
}

/* End of code generation (SmRG_mixtureModelFitting.c) */
