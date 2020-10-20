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
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_newPost.h"
#include "SmRG_mixtureModelFitting_newPost_emxutil.h"
#include "combineVectorElements.h"
#include "mean.h"
#include "sum.h"
#include "psi.h"
#include "exp.h"
#include "gammaln.h"
#include "SmRG_normlike.h"

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
  *K0_double, double *vB_double, double *mu_sk, double *rk, emxArray_real_T
  *p_tot, emxArray_real_T *a1)
{
  int i0;
  int loop_ub;
  emxArray_uint32_T *i_ok;
  unsigned int a_idx_0;
  double pk;
  double pb;
  double pa;
  int varargin_1;
  int n_of_pixel;
  double log_lik_a[5000];
  double log_lik_b[5000];
  int count;
  double tol_nbin;
  double tol_gauss;
  emxArray_real_T *psi_a;
  emxArray_real_T *psi_b;
  emxArray_real_T *dk;
  emxArray_real_T *b;
  emxArray_real_T *a;
  emxArray_real_T *b_b;
  emxArray_real_T *b_a1;
  emxArray_real_T *c_a1;
  boolean_T exitg1;
  emxArray_int32_T *r0;
  double b_K0_double[2];
  double psi_b_tmp;
  double B_tmp_data[1];
  int B_tmp_size[2];
  int nx;
  double d0;
  double d1;
  double d2;
  emxArray_real_T d_a1;
  int c_b[1];

  /*  SmRG_mixtureModelFitting: */
  /*            fits model described in [1] on data input Vin. The model */
  /*            describes a single class k of signal pixels (negative-binomial) */
  /*            and the background (gaussian).For more details refer to [1]. */
  /*            For negative binomial parameter update refer to [2]. */
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
    a_idx_0 = (unsigned int)(Vin->size[0] * Vin->size[1]);
    i0 = i_ok->size[0] * i_ok->size[1];
    i_ok->size[0] = 1;
    i_ok->size[1] = (int)a_idx_0;
    emxEnsureCapacity_uint32_T(i_ok, i0);
    loop_ub = (int)a_idx_0 - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      i_ok->data[i0] = 1U + i0;
    }
  }

  pk = *mu_sk / ((*mu_sk * *mu_sk + *mu_sk * *rk) / *rk);

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
  emxInit_real_T(&b_b, 1);
  emxInit_real_T(&b_a1, 1);
  emxInit_real_T(&c_a1, 2);
  exitg1 = false;
  while ((!exitg1) && ((tol_gauss > 0.0001) && (tol_nbin > 0.0001))) {
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
      /*  K0_double = double(K0); vB_double = double(vB); */
      /*   Gaussian log likelihood */
      b_K0_double[0] = *K0_double;
      b_K0_double[1] = sqrt(*vB_double);
      psi_b_tmp = -SmRG_normlike(b_K0_double, Vin->data[loop_ub]);

      /*  USEFUL IN DEBUG MODE */
      /*  if isnan(psi_b_tmp);keyboard;end */
      /* if (Vtmp-K0_double)<0;K0_double=double(min(V));end */
      psi_b->data[loop_ub] = psi_b_tmp;

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
      /*  if min(n,m) > 1 */
      /*      error(message('stats:nbinlike:InvalidData')); */
      /*  end */
      tol_nbin = 1.0 / (1.0 / *rk);
      tol_gauss = tol_nbin / (tol_nbin + *mu_sk);

      /*  Out of range or missing parameters or data return NaN.  Infinite */
      /*  values for R correspond to a Poisson, but its mean cannot be determined */
      /*  from the (R,P) parametrization. */
      if ((!(0.0 < tol_nbin)) || (rtIsInf(tol_nbin) || rtIsNaN(tol_nbin)) ||
          (!(0.0 < tol_gauss)) || (!(tol_gauss <= 1.0)) || ((!(0.0 <= Vin->
             data[loop_ub])) || (rtIsInf(Vin->data[loop_ub]) || rtIsNaN
            (Vin->data[loop_ub])) || (!(Vin->data[loop_ub] == rt_roundd_snf
             (Vin->data[loop_ub]))))) {
        tol_nbin = rtNaN;
      } else {
        d0 = tol_nbin + Vin->data[loop_ub];
        gammaln(&d0);
        d1 = Vin->data[loop_ub] + 1.0;
        gammaln(&d1);
        d2 = tol_nbin;
        gammaln(&d2);
        tol_nbin = -(((d0 - d1) - d2) + tol_nbin * log(tol_gauss)) - Vin->
          data[loop_ub] * log(1.0 - tol_gauss);
      }

      /*  USEFUL IN DEBUG MODE */
      /*   if isnan(psi_a_tmp);keyboard;end */
      psi_a->data[loop_ub] = -tol_nbin;

      /*  find who's bigger between [log(P(x|a))+log(P(a)] e [log(P(x|b))+log(P(b)] */
      /*  to use the logarithm property of sum: log(sum(ai))=log_b(a0)+log_b(1+sum(b^(log_b(ai)-log_b(a0)) */
      /*  with a0>a1>...>aN */
      /*  in this way it is possible to evaluate the log posterior without */
      /*  calculating the exponential: this helps for numerical stability */
      /*  (sometimes exponential blows up) */
      /*  find the denominator of log posterior */
      if (-tol_nbin + log(pa) > psi_b_tmp + log(pb)) {
        tol_gauss = (-tol_nbin + log(pa)) + log(1.0 + exp(((psi_b_tmp + log(pb))
          - (-tol_nbin)) + log(pa)));
      } else {
        tol_gauss = (psi_b_tmp + log(pb)) + log(1.0 + exp(((-tol_nbin + log(pa))
          - psi_b_tmp) + log(pb)));
      }

      /*  old way to evaluate the log posterior */
      /*   b(ii,1) = ((psi_b_tmp)+log(pb))-log(exp(psi_b_tmp)*pb+exp(psi_a_tmp)*pa); */
      /*   a(ii,1) = ((psi_a_tmp)+log(pa))-log(exp(psi_b_tmp)*pb+exp(psi_a_tmp)*pa); */
      /*  new way; actually much better */
      b->data[loop_ub] = (psi_b_tmp + log(pb)) - tol_gauss;
      a->data[loop_ub] = (-tol_nbin + log(pa)) - tol_gauss;
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

    /*  USEFUL IN DEBUG MODE */
    /*   if sum(isnan(b1))~=0;keyboard;end */
    /*   if sum(isnan(a1))~=0;keyboard;end */
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
    i0 = b_b->size[0];
    b_b->size[0] = Vin->size[0] * Vin->size[1];
    emxEnsureCapacity_real_T(b_b, i0);
    loop_ub = Vin->size[0] * Vin->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_b->data[i0] = Vin->data[i0] - *K0_double;
    }

    a_idx_0 = (unsigned int)b_b->size[0];
    i0 = dk->size[0];
    dk->size[0] = (int)a_idx_0;
    emxEnsureCapacity_real_T(dk, i0);
    a_idx_0 = (unsigned int)b_b->size[0];
    nx = (int)a_idx_0;
    for (loop_ub = 0; loop_ub < nx; loop_ub++) {
      dk->data[loop_ub] = b_b->data[loop_ub] * b_b->data[loop_ub];
    }

    loop_ub = b->size[0];
    i0 = dk->size[0];
    dk->size[0] = loop_ub;
    emxEnsureCapacity_real_T(dk, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      dk->data[i0] *= b->data[i0];
    }

    *vB_double = sum(dk) / B_tmp_data[0];

    /*  update Neg-Bin Parameters according to [2] */
    i0 = b_b->size[0];
    b_b->size[0] = Vin->size[0] * Vin->size[1];
    emxEnsureCapacity_real_T(b_b, i0);
    loop_ub = Vin->size[0] * Vin->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_b->data[i0] = *rk + Vin->data[i0];
    }

    nx = b_b->size[0];
    a_idx_0 = (unsigned int)b_b->size[0];
    i0 = dk->size[0];
    dk->size[0] = (int)a_idx_0;
    emxEnsureCapacity_real_T(dk, i0);
    for (loop_ub = 0; loop_ub < nx; loop_ub++) {
      dk->data[loop_ub] = minusDoubleScalarPsi(b_b->data[loop_ub]);
    }

    tol_nbin = minusDoubleScalarPsi(*rk);
    i0 = dk->size[0];
    emxEnsureCapacity_real_T(dk, i0);
    loop_ub = dk->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      dk->data[i0] = *rk * (dk->data[i0] - tol_nbin);
    }

    tol_nbin = (1.0 - 1.0 / (1.0 - pk)) - 1.0 / log(pk);
    loop_ub = a1->size[0];
    i0 = b_b->size[0];
    b_b->size[0] = loop_ub;
    emxEnsureCapacity_real_T(b_b, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_b->data[i0] = a1->data[i0] * dk->data[i0];
    }

    loop_ub = a1->size[0];
    i0 = b_a1->size[0];
    b_a1->size[0] = loop_ub;
    emxEnsureCapacity_real_T(b_a1, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_a1->data[i0] = a1->data[i0] * (Vin->data[i0] - (1.0 - tol_nbin) *
        dk->data[i0]);
    }

    pk = tol_nbin * sum(b_b) / sum(b_a1);
    loop_ub = a1->size[0];
    i0 = dk->size[0];
    dk->size[0] = loop_ub;
    emxEnsureCapacity_real_T(dk, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      dk->data[i0] *= a1->data[i0];
    }

    loop_ub = a1->size[0];
    d_a1 = *a1;
    c_b[0] = loop_ub;
    d_a1.size = &c_b[0];
    d_a1.numDimensions = 1;
    *rk = -(sum(dk) / sum(&d_a1)) / log(pk);

    /*  USEFUL IN DEBUG MODE */
    /*   if ~isreal(pk) */
    /*       keyboard */
    /*   end */
    *mu_sk = pk * (*rk * (1.0 - pk) / (pk * pk));

    /*  today's posterior becomes tomorrow's prior :) */
    loop_ub = a1->size[0];
    i0 = c_a1->size[0] * c_a1->size[1];
    c_a1->size[0] = 1;
    c_a1->size[1] = loop_ub;
    emxEnsureCapacity_real_T(c_a1, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_a1->data[i0] = a1->data[i0];
    }

    pa = mean(c_a1);
    loop_ub = b->size[0];
    i0 = c_a1->size[0] * c_a1->size[1];
    c_a1->size[0] = 1;
    c_a1->size[1] = loop_ub;
    emxEnsureCapacity_real_T(c_a1, i0);
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_a1->data[i0] = b->data[i0];
    }

    pb = mean(c_a1);

    /*  update for tollerance check */
    count++;
    log_lik_a[count] = combineVectorElements(psi_a) / (double)psi_a->size[0];
    log_lik_b[count] = combineVectorElements(psi_b) / (double)psi_b->size[0];

    /*  prevent infinite loop */
    if (count + 1 > 4999) {
      exitg1 = true;
    } else {
      /*  reset tol_nbin and tol_gauss */
      tol_nbin = fabs(log_lik_a[count] - log_lik_a[count - 1]) / fabs
        (log_lik_a[count]);

      /*  update for tollerance check gaussian */
      log_lik_b[count] = combineVectorElements(psi_b) / (double)psi_b->size[0];
      tol_gauss = fabs(log_lik_b[count] - log_lik_b[count - 1]) / fabs
        (log_lik_b[count]);
    }
  }

  emxFree_real_T(&c_a1);
  emxFree_real_T(&b_a1);
  emxFree_real_T(&b_b);
  emxFree_real_T(&a);
  emxFree_real_T(&b);
  emxFree_real_T(&dk);
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

/* End of code generation (SmRG_mixtureModelFitting_newPost.c) */
