/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_mixtureModelFitting_newPost.c
 *
 * Code generation for function 'SmRG_mixtureModelFitting_newPost'
 *
 */

/* Include files */
#include "SmRG_mixtureModelFitting_newPost.h"
#include "SmRG_mixtureModelFitting_newPost_data.h"
#include "SmRG_mixtureModelFitting_newPost_emxutil.h"
#include "SmRG_mixtureModelFitting_newPost_initialize.h"
#include "SmRG_mixtureModelFitting_newPost_types.h"
#include "gammaln.h"
#include "psi.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

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

void SmRG_mixtureModelFitting_newPost(const emxArray_real_T *Vin, double
  *K0_double, double *vB_double, double *mu_sk, double *rk, double tol,
  emxArray_real_T *p_tot, emxArray_real_T *a1)
{
  emxArray_int32_T *r;
  emxArray_real_T *a;
  emxArray_real_T *b;
  emxArray_real_T *dk;
  emxArray_real_T *psi_a;
  emxArray_real_T *psi_b;
  emxArray_real_T *x;
  emxArray_uint32_T *i_ok;
  double log_lik_a[5000];
  double log_lik_b[5000];
  double d;
  double d1;
  double d2;
  double d3;
  double nlogL;
  double pa;
  double pb;
  double pk;
  double tol_gauss;
  double tol_nbin;
  double y;
  int count;
  int k;
  int n_of_pixel;
  int nx;
  boolean_T exitg1;
  if (!isInitialized_SmRG_mixtureModelFitting_newPost) {
    SmRG_mixtureModelFitting_newPost_initialize();
  }

  /*  SmRG_mixtureModelFitting: */
  /*            fits model described in [1] on data input Vin. The model */
  /*            describes a single class k of signal pixels (negative-binomial) */
  /*            and the background (gaussian).For more details refer to [1]. */
  /*            For negative binomial parameter update refer to [2]. */
  /*  */
  /*  Syntax: */
  /*            [p_tot,a1]=SmRG_mixtureModelFitting_newPost(vin) */
  /*            [p_tot,a1,K0_double,vB,mu_sk,rk]=SmRG_mixtureModelFitting_newPost(vin) */
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
  /* [2]: Chunmao Huang et al 2019 J. Phys.: Conf. Ser. 1324 012093 */
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
  k = p_tot->size[0] * p_tot->size[1];
  p_tot->size[0] = 1;
  p_tot->size[1] = Vin->size[1];
  emxEnsureCapacity_real_T(p_tot, k);
  nx = Vin->size[1];
  for (k = 0; k < nx; k++) {
    p_tot->data[k] = 0.0;
  }

  k = a1->size[0] * a1->size[1];
  a1->size[0] = 1;
  a1->size[1] = Vin->size[1];
  emxEnsureCapacity_real_T(a1, k);
  nx = Vin->size[1];
  for (k = 0; k < nx; k++) {
    a1->data[k] = 0.0;
  }

  emxInit_uint32_T(&i_ok, 2);
  if (Vin->size[1] < 1) {
    i_ok->size[0] = 1;
    i_ok->size[1] = 0;
  } else {
    k = i_ok->size[0] * i_ok->size[1];
    i_ok->size[0] = 1;
    i_ok->size[1] = Vin->size[1];
    emxEnsureCapacity_uint32_T(i_ok, k);
    nx = Vin->size[1] - 1;
    for (k = 0; k <= nx; k++) {
      i_ok->data[k] = k + 1U;
    }
  }

  pk = *mu_sk / ((*mu_sk * *mu_sk + *mu_sk * *rk) / *rk);

  /*  Mixing prior probability (0.5 0.5) */
  pb = 0.5;
  pa = 0.5;

  /*  prior */
  /*  EM initial conditions */
  n_of_pixel = Vin->size[1];
  memset(&log_lik_a[0], 0, 5000U * sizeof(double));
  memset(&log_lik_b[0], 0, 5000U * sizeof(double));
  log_lik_a[0] = 0.0;
  log_lik_b[0] = 0.0;
  log_lik_a[1] = 10.0;
  log_lik_b[1] = 10.0;

  /*  tol   = 1e-4; */
  count = 1;

  /*  reset tol_nbin and tol_gauss */
  tol_nbin = 1.0;

  /*  update for tollerance check gaussian */
  tol_gauss = 1.0;

  /*  USEFUL IN DEBUG MODE: PLOTS THE MODEL FITTING ON DATA HIST */
  /*      figure */
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
  emxInit_real_T(&dk, 1);
  emxInit_real_T(&b, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&x, 1);
  exitg1 = false;
  while ((!exitg1) && ((tol_gauss > tol) && (tol_nbin > tol))) {
    /*  initialize log likelihood */
    k = psi_a->size[0];
    psi_a->size[0] = n_of_pixel;
    emxEnsureCapacity_real_T(psi_a, k);
    k = psi_b->size[0];
    psi_b->size[0] = n_of_pixel;
    emxEnsureCapacity_real_T(psi_b, k);

    /*  initialize log posterior */
    k = b->size[0] * b->size[1];
    b->size[0] = n_of_pixel;
    b->size[1] = 1;
    emxEnsureCapacity_real_T(b, k);
    k = a->size[0] * a->size[1];
    a->size[0] = n_of_pixel;
    a->size[1] = 1;
    emxEnsureCapacity_real_T(a, k);
    for (k = 0; k < n_of_pixel; k++) {
      psi_a->data[k] = 0.0;
      psi_b->data[k] = 0.0;
      b->data[k] = 0.0;
      a->data[k] = 0.0;
    }

    /*  E-step */
    if (0 <= n_of_pixel - 1) {
      d = sqrt(*vB_double);
    }

    for (nx = 0; nx < n_of_pixel; nx++) {
      /*  cast */
      /*  K0_double = double(K0); vB_double = double(vB); */
      /*   Gaussian log likelihood */
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
      tol_nbin = d;
      if (d <= 0.0) {
        tol_nbin = rtNaN;
      }

      /*  Compute the individual log-likelihood terms. */
      d1 = Vin->data[nx];
      tol_gauss = (d1 - *K0_double) / tol_nbin;

      /*  Sum up the individual contributions, and return the negative log-likelihood. */
      nlogL = -(-0.5 * tol_gauss * tol_gauss - log(2.5066282746310002 * tol_nbin));

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
      /*  USEFUL IN DEBUG MODE */
      /*  if isnan(psi_b_tmp);keyboard;end */
      /* if (Vtmp-K0_double)<0;K0_double=double(min(V));end */
      psi_b->data[nx] = -nlogL;

      /*  USEFUL IN DEBUG MODE */
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
      /*  if nargin < 2,  */
      /*      error(message('stats:nbinlike:TooFewInputs'));  */
      /*  end */
      /*   */
      /*  if min(n,m) > 1 */
      /*      error(message('stats:nbinlike:InvalidData')); */
      /*  end */
      /*   */
      /*  if nargout == 2 & max(m,n) == 1 */
      /*      error(message('stats:nbinlike:NotEnoughData')); */
      /*  end */
      tol_gauss = 1.0 / (1.0 / *rk);
      tol_nbin = tol_gauss / (tol_gauss + *mu_sk);

      /*  Out of range or missing parameters or data return NaN.  Infinite */
      /*  values for R correspond to a Poisson, but its mean cannot be determined */
      /*  from the (R,P) parametrization. */
      if ((!(0.0 < tol_gauss)) || (rtIsInf(tol_gauss) || rtIsNaN(tol_gauss)) ||
          (!(0.0 < tol_nbin)) || (!(tol_nbin <= 1.0)) || ((!(0.0 <= d1)) ||
           (rtIsInf(d1) || rtIsNaN(d1)) || (!(d1 == rt_roundd_snf(d1))))) {
        tol_gauss = rtNaN;
      } else {
        y = tol_gauss + d1;
        gammaln(&y);
        d2 = d1 + 1.0;
        gammaln(&d2);
        d3 = tol_gauss;
        gammaln(&d3);
        tol_gauss = -(((y - d2) - d3) + tol_gauss * log(tol_nbin)) - d1 * log
          (1.0 - tol_nbin);
      }

      /*  USEFUL IN DEBUG MODE */
      /*   if isnan(psi_a_tmp);keyboard;end */
      psi_a->data[nx] = -tol_gauss;

      /*  find who's bigger between [log(P(x|a))+log(P(a)] e [log(P(x|b))+log(P(b)] */
      /*  to use the logarithm property of sum: log(sum(ai))=log_b(a0)+log_b(1+sum(b^(log_b(ai)-log_b(a0)) */
      /*  with a0>a1>...>aN */
      /*  in this way it is possible to evaluate the log posterior without */
      /*  calculating the exponential: this helps for numerical stability */
      /*  (sometimes exponential blows up) */
      /*  find the denominator of log posterior */
      if (-tol_gauss + log(pa) > -nlogL + log(pb)) {
        tol_nbin = (-tol_gauss + log(pa)) + log(exp(((-nlogL + log(pb)) -
          (-tol_gauss)) + log(pa)) + 1.0);
      } else {
        tol_nbin = (-nlogL + log(pb)) + log(exp(((-tol_gauss + log(pa)) -
          (-nlogL)) + log(pb)) + 1.0);
      }

      /*  old way to evaluate the log posterior */
      /*   b(ii,1) = ((psi_b_tmp)+log(pb))-log(exp(psi_b_tmp)*pb+exp(psi_a_tmp)*pa); */
      /*   a(ii,1) = ((psi_a_tmp)+log(pa))-log(exp(psi_b_tmp)*pb+exp(psi_a_tmp)*pa); */
      /*  new way; actually much better */
      b->data[nx] = (-nlogL + log(pb)) - tol_nbin;
      a->data[nx] = (-tol_gauss + log(pa)) - tol_nbin;
    }

    nx = b->size[0];
    for (k = 0; k < nx; k++) {
      b->data[k] = exp(b->data[k]);
    }

    nx = a->size[0];
    for (k = 0; k < nx; k++) {
      a->data[k] = exp(a->data[k]);
    }

    k = a1->size[0] * a1->size[1];
    a1->size[0] = a->size[0];
    a1->size[1] = 1;
    emxEnsureCapacity_real_T(a1, k);
    nx = a->size[0] * a->size[1];
    for (k = 0; k < nx; k++) {
      a1->data[k] = a->data[k];
    }

    /*  USEFUL IN DEBUG MODE */
    /*   if sum(isnan(b1))~=0;keyboard;end */
    /*   if sum(isnan(a1))~=0;keyboard;end */
    /*  M-step */
    /*  update Gaussian Parameters */
    nx = b->size[0];
    if (b->size[0] == 0) {
      tol_gauss = 0.0;
    } else {
      tol_gauss = b->data[0];
      for (k = 2; k <= nx; k++) {
        tol_gauss += b->data[k - 1];
      }
    }

    k = x->size[0];
    x->size[0] = b->size[0];
    emxEnsureCapacity_real_T(x, k);
    nx = b->size[0];
    for (k = 0; k < nx; k++) {
      x->data[k] = b->data[k] * Vin->data[k];
    }

    nx = x->size[0];
    if (x->size[0] == 0) {
      y = 0.0;
    } else {
      y = x->data[0];
      for (k = 2; k <= nx; k++) {
        y += x->data[k - 1];
      }
    }

    *K0_double = y / tol_gauss;
    k = dk->size[0];
    dk->size[0] = Vin->size[1];
    emxEnsureCapacity_real_T(dk, k);
    nx = Vin->size[1];
    for (k = 0; k < nx; k++) {
      dk->data[k] = Vin->data[k] - *K0_double;
    }

    k = x->size[0];
    x->size[0] = dk->size[0];
    emxEnsureCapacity_real_T(x, k);
    nx = dk->size[0];
    for (k = 0; k < nx; k++) {
      x->data[k] = dk->data[k] * dk->data[k];
    }

    nx = b->size[0];
    k = x->size[0];
    x->size[0] = b->size[0];
    emxEnsureCapacity_real_T(x, k);
    for (k = 0; k < nx; k++) {
      x->data[k] *= b->data[k];
    }

    nx = x->size[0];
    if (x->size[0] == 0) {
      y = 0.0;
    } else {
      y = x->data[0];
      for (k = 2; k <= nx; k++) {
        y += x->data[k - 1];
      }
    }

    *vB_double = y / tol_gauss;

    /*  update Neg-Bin Parameters according to [2] */
    k = x->size[0];
    x->size[0] = Vin->size[1];
    emxEnsureCapacity_real_T(x, k);
    nx = Vin->size[1];
    for (k = 0; k < nx; k++) {
      x->data[k] = *rk + Vin->data[k];
    }

    nx = x->size[0];
    k = dk->size[0];
    dk->size[0] = x->size[0];
    emxEnsureCapacity_real_T(dk, k);
    for (k = 0; k < nx; k++) {
      dk->data[k] = minusDoubleScalarPsi(x->data[k]);
    }

    y = minusDoubleScalarPsi(*rk);
    nx = dk->size[0];
    for (k = 0; k < nx; k++) {
      dk->data[k] = *rk * (dk->data[k] - y);
    }

    k = x->size[0];
    x->size[0] = a->size[0];
    emxEnsureCapacity_real_T(x, k);
    nx = a->size[0];
    for (k = 0; k < nx; k++) {
      x->data[k] = a->data[k] * dk->data[k];
    }

    nx = x->size[0];
    if (x->size[0] == 0) {
      tol_nbin = 0.0;
    } else {
      tol_nbin = x->data[0];
      for (k = 2; k <= nx; k++) {
        tol_nbin += x->data[k - 1];
      }
    }

    tol_gauss = (1.0 - 1.0 / (1.0 - pk)) - 1.0 / log(pk);
    nx = a->size[0];
    k = dk->size[0];
    dk->size[0] = a->size[0];
    emxEnsureCapacity_real_T(dk, k);
    for (k = 0; k < nx; k++) {
      dk->data[k] = a->data[k] * (Vin->data[k] - (1.0 - tol_gauss) * dk->data[k]);
    }

    nx = dk->size[0];
    if (dk->size[0] == 0) {
      y = 0.0;
    } else {
      y = dk->data[0];
      for (k = 2; k <= nx; k++) {
        y += dk->data[k - 1];
      }
    }

    pk = tol_gauss * tol_nbin / y;
    nx = a->size[0];
    if (a->size[0] == 0) {
      y = 0.0;
    } else {
      y = a->data[0];
      for (k = 2; k <= nx; k++) {
        y += a->data[k - 1];
      }
    }

    *rk = -(tol_nbin / y) / log(pk);

    /*  USEFUL IN DEBUG MODE */
    /*   if ~isreal(pk) */
    /*       keyboard */
    /*   end */
    *mu_sk = pk * (*rk * (1.0 - pk) / (pk * pk));

    /*  today's posterior becomes tomorrow's prior :) */
    nx = a->size[0];
    if (a->size[0] == 0) {
      y = 0.0;
    } else {
      y = a->data[0];
      for (k = 2; k <= nx; k++) {
        y += a->data[k - 1];
      }
    }

    pa = y / (double)a->size[0];
    nx = b->size[0];
    if (b->size[0] == 0) {
      y = 0.0;
    } else {
      y = b->data[0];
      for (k = 2; k <= nx; k++) {
        y += b->data[k - 1];
      }
    }

    pb = y / (double)b->size[0];

    /*  update for tollerance check */
    count++;
    nx = psi_a->size[0];
    if (psi_a->size[0] == 0) {
      y = 0.0;
    } else {
      y = psi_a->data[0];
      for (k = 2; k <= nx; k++) {
        y += psi_a->data[k - 1];
      }
    }

    log_lik_a[count] = y / (double)psi_a->size[0];
    nx = psi_b->size[0];
    if (psi_b->size[0] == 0) {
      y = 0.0;
    } else {
      y = psi_b->data[0];
      for (k = 2; k <= nx; k++) {
        y += psi_b->data[k - 1];
      }
    }

    tol_gauss = y / (double)psi_b->size[0];
    log_lik_b[count] = tol_gauss;

    /*  prevent infinite loop */
    if (count + 1 > 4999) {
      exitg1 = true;
    } else {
      /*  reset tol_nbin and tol_gauss */
      tol_nbin = fabs(log_lik_a[count] - log_lik_a[count - 1]) / fabs
        (log_lik_a[count]);

      /*  update for tollerance check gaussian */
      log_lik_b[count] = tol_gauss;
      tol_gauss = fabs(log_lik_b[count] - log_lik_b[count - 1]) / fabs
        (log_lik_b[count]);
    }
  }

  emxFree_real_T(&x);
  emxFree_real_T(&a);
  emxFree_real_T(&b);
  emxFree_real_T(&dk);
  emxFree_real_T(&psi_b);
  emxFree_real_T(&psi_a);
  emxInit_int32_T(&r, 2);
  k = r->size[0] * r->size[1];
  r->size[0] = 1;
  r->size[1] = i_ok->size[1];
  emxEnsureCapacity_int32_T(r, k);
  nx = i_ok->size[0] * i_ok->size[1];
  for (k = 0; k < nx; k++) {
    r->data[k] = (int)i_ok->data[k];
  }

  emxFree_uint32_T(&i_ok);
  nx = r->size[0] * r->size[1];
  for (k = 0; k < nx; k++) {
    p_tot->data[r->data[k] - 1] = a1->data[k];
  }

  emxFree_int32_T(&r);

  /*  p_tot(i_nan)=1; */
}

/* End of code generation (SmRG_mixtureModelFitting_newPost.c) */
