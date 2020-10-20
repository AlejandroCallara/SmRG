/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_mixtureModelFitting_multmix.c
 *
 * Code generation for function 'SmRG_mixtureModelFitting_multmix'
 *
 */

/* Include files */
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "SmRG_mixtureModelFitting_multmix_emxutil.h"
#include "rdivide_helper.h"
#include "abs.h"
#include "power.h"
#include "sum.h"
#include "mean.h"
#include "structConstructorHelper.h"
#include "repmat.h"
#include "SmRG_kmeans_opt.h"
#include "log.h"
#include "psi.h"
#include "exp.h"
#include "sort1.h"
#include "nbinlike_mu.h"
#include "SmRG_normlike.h"
#include "sqrt.h"
#include <stdio.h>

/* Type Definitions */
#ifndef struct_c_emxArray_sQOWHQUVurFQwFR6JFyF
#define struct_c_emxArray_sQOWHQUVurFQwFR6JFyF

struct c_emxArray_sQOWHQUVurFQwFR6JFyF
{
  cell_wrap_5 data[4];
  int size[1];
};

#endif                                 /*struct_c_emxArray_sQOWHQUVurFQwFR6JFyF*/

#ifndef typedef_emxArray_cell_wrap_5_4
#define typedef_emxArray_cell_wrap_5_4

typedef struct c_emxArray_sQOWHQUVurFQwFR6JFyF emxArray_cell_wrap_5_4;

#endif                                 /*typedef_emxArray_cell_wrap_5_4*/

#ifndef struct_c_emxArray_spWy85dGSaZ7aeYViKME
#define struct_c_emxArray_spWy85dGSaZ7aeYViKME

struct c_emxArray_spWy85dGSaZ7aeYViKME
{
  struct_T data[4];
  int size[1];
};

#endif                                 /*struct_c_emxArray_spWy85dGSaZ7aeYViKME*/

#ifndef typedef_emxArray_struct_T_4
#define typedef_emxArray_struct_T_4

typedef struct c_emxArray_spWy85dGSaZ7aeYViKME emxArray_struct_T_4;

#endif                                 /*typedef_emxArray_struct_T_4*/

/* Function Definitions */
void SmRG_mixtureModelFitting_multmix(const emxArray_real_T *Vin, double
  *K0_double, double *vB_double, double mu_sk_data[], int mu_sk_size[2], double
  rk_data[], int rk_size[2], emxArray_real_T *p_tot, emxArray_real_T *a1)
{
  emxArray_real_T *a;
  double nnegbin;
  int varargin_2;
  int i7;
  int loop_ub;
  emxArray_real_T *b;
  double tol_gauss;
  double mu_sk_tmp[4];
  double rk_tmp[4];
  emxArray_real_T *b_psort;
  emxArray_real_T b_mu_sk_tmp;
  int iv0[2];
  int psort_size[2];
  double psort_data[4];
  double v_sk_data[4];
  int v_sk_size[2];
  static int iv1[2] = { 1, 4 };

  double pk_data[4];
  int pk_size[2];
  double pb;
  int meana1_size[2];
  double meana1_data[4];
  emxArray_real_T *log_lik_a;
  int b_varargin_2;
  double log_lik_b[5000];
  int count;
  double tmp_data[4];
  int tmp_size[2];
  double b_v_sk_data[4];
  int b_v_sk_size[2];
  int n;
  int idx;
  double tol_nbin;
  int k;
  boolean_T exitg1;
  double d4;
  emxArray_real_T *b1;
  emxArray_real_T *psi_a;
  emxArray_real_T *psi_b;
  emxArray_real_T *dk;
  emxArray_real_T *b_b1;
  emxArray_real_T *b_a1;
  emxArray_real_T *c_a1;
  emxArray_real_T *d_a1;
  emxArray_int32_T *c_psort;
  int exitg2;
  int loop_ub_tmp;
  double b_K0_double[2];
  double psi_b_tmp;
  int b_mu_sk_size[2];
  double b_mu_sk_data[8];
  int v_opt_size[1];
  double C_data[8];
  int SUMD_size[1];
  cell_wrap_5 r0;
  cell_wrap_5 rv0[1];
  emxArray_cell_wrap_5_4 r1;
  int i8;
  emxArray_struct_T_4 b_idx;
  int mus_tmp_size[2];
  int mus_tmp_size_tmp;
  double mus_tmp_data[3];
  double sigmas_tmp_data[3];
  double p_tmp_data[3];
  signed char b_tmp_data[4];
  emxArray_real_T b_mus_tmp_data;
  double b_sigmas_tmp_data[3];
  signed char c_tmp_data[4];
  int rk[1];
  int log_lik_a_size[2];
  signed char d_tmp_data[4];
  emxArray_real_T c_v_sk_data;
  emxArray_real_T b_pk_data;
  boolean_T exitg3;
  emxInit_real_T(&a, 2);

  /*            fits model described in [1] on data input Vin. The model */
  /*            describes a single class k of signal pixels (negative-binomial) */
  /*            and the background (gaussian).For more details refer to [1]. */
  /*            For negative binomial parameter update refer to [2]. */
  /*  */
  /*  Syntax: */
  /*            [p_tot,a1]=SmRG_mixtureModelFitting_multimix(vin) */
  /*            [p_tot,a1,K0_double,vB,mu_sk,rk]=SmRG_mixtureModelFitting_multimix(vin) */
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
  /*  */
  /*  % if nargin<2 */
  /*  %     % if background not specified never do check on background */
  /*  %     background = 0; */
  /*  % end */
  /*  */
  /*  if sum(Vin(:)<0) */
  /*      error ('only images with positive values are allowed') */
  /*  end */
  /*  if ~isrow(Vin) */
  /*      Vin = Vin'; */
  /*  end */
  /*  if ~isrow(Vin) */
  /*      mu_sk = mu_sk'; */
  /*  end,mu_sk,rk */
  nnegbin = 4.0;

  /*  check for nans in data */
  varargin_2 = Vin->size[1];

  /*  pre-allocate p_tot */
  /*  p_tot=zeros(1,l); */
  i7 = a->size[0] * a->size[1];
  a->size[0] = Vin->size[1];
  a->size[1] = 4;
  emxEnsureCapacity_real_T(a, i7);
  loop_ub = Vin->size[1] << 2;
  for (i7 = 0; i7 < loop_ub; i7++) {
    a->data[i7] = 0.0;
  }

  i7 = a1->size[0] * a1->size[1];
  a1->size[0] = Vin->size[1];
  a1->size[1] = 4;
  emxEnsureCapacity_real_T(a1, i7);
  loop_ub = Vin->size[1] << 2;
  for (i7 = 0; i7 < loop_ub; i7++) {
    a1->data[i7] = 0.0;
  }

  emxInit_real_T(&b, 1);
  i7 = b->size[0];
  b->size[0] = Vin->size[1];
  emxEnsureCapacity_real_T(b, i7);
  loop_ub = Vin->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    b->data[i7] = 0.0;
  }

  /*  if no initial conditions then do it data driven */
  /* else get them */
  tol_gauss = mean(Vin);
  mu_sk_tmp[0] = tol_gauss;
  mu_sk_tmp[1] = tol_gauss;
  mu_sk_tmp[2] = tol_gauss;
  mu_sk_tmp[3] = tol_gauss;

  /*  mean */
  if (1 > mu_sk_size[1]) {
    loop_ub = 0;
  } else {
    loop_ub = mu_sk_size[1];
  }

  if (0 <= loop_ub - 1) {
    memcpy(&mu_sk_tmp[0], &mu_sk_data[0], (unsigned int)(loop_ub * (int)sizeof
            (double)));
  }

  rk_tmp[0] = 10.0;
  rk_tmp[1] = 10.0;
  rk_tmp[2] = 10.0;
  rk_tmp[3] = 10.0;
  if (1 > mu_sk_size[1]) {
    loop_ub = 0;
  } else {
    loop_ub = mu_sk_size[1];
  }

  if (0 <= loop_ub - 1) {
    memcpy(&rk_tmp[0], &rk_data[0], (unsigned int)(loop_ub * (int)sizeof(double)));
  }

  mu_sk_size[0] = 1;
  mu_sk_size[1] = 4;
  mu_sk_data[0] = mu_sk_tmp[0];
  mu_sk_data[1] = mu_sk_tmp[1];
  mu_sk_data[2] = mu_sk_tmp[2];
  mu_sk_data[3] = mu_sk_tmp[3];
  rk_size[0] = 1;
  rk_size[1] = 4;
  rk_data[0] = rk_tmp[0];
  rk_data[1] = rk_tmp[1];
  rk_data[2] = rk_tmp[2];
  rk_data[3] = rk_tmp[3];
  emxInit_real_T(&b_psort, 2);
  b_mu_sk_tmp.numDimensions = 2;
  iv0[0] = 1;
  iv0[1] = 4;
  b_mu_sk_tmp.size = &iv0[0];
  b_mu_sk_tmp.allocatedSize = 4;
  b_mu_sk_tmp.data = &mu_sk_tmp[0];
  power(&b_mu_sk_tmp, b_psort);
  psort_size[0] = 1;
  psort_size[1] = b_psort->size[1];
  loop_ub = b_psort->size[0] * b_psort->size[1];
  for (i7 = 0; i7 < loop_ub; i7++) {
    psort_data[i7] = b_psort->data[i7] + mu_sk_tmp[i7] * rk_tmp[i7];
  }

  rdivide_helper(psort_data, psort_size, rk_tmp, v_sk_data, v_sk_size);
  rdivide_helper(mu_sk_tmp, iv1, v_sk_data, pk_data, pk_size);

  /*  Mixing prior probability (0.5 0.5/nnegbin) */
  pb = 0.5;
  meana1_size[1] = 4;
  meana1_data[0] = 0.125;
  meana1_data[1] = 0.125;
  meana1_data[2] = 0.125;
  meana1_data[3] = 0.125;
  emxInit_real_T(&log_lik_a, 2);

  /*  EM initial conditions */
  b_varargin_2 = Vin->size[1];
  i7 = log_lik_a->size[0] * log_lik_a->size[1];
  log_lik_a->size[0] = 5000;
  log_lik_a->size[1] = 4;
  emxEnsureCapacity_real_T(log_lik_a, i7);
  for (i7 = 0; i7 < 20000; i7++) {
    log_lik_a->data[i7] = 0.0;
  }

  memset(&log_lik_b[0], 0, 5000U * sizeof(double));
  log_lik_b[0] = 0.0;
  log_lik_b[1] = 10.0;
  log_lik_a->data[0] = 0.0;
  log_lik_a->data[1] = 10.0;
  log_lik_a->data[5000] = 0.0;
  log_lik_a->data[5001] = 10.0;
  log_lik_a->data[10000] = 0.0;
  log_lik_a->data[10001] = 10.0;
  log_lik_a->data[15000] = 0.0;
  log_lik_a->data[15001] = 10.0;
  count = 1;

  /*  reset tol_nbin and tol_gauss */
  psort_size[0] = 1;
  psort_size[1] = 4;
  psort_data[0] = log_lik_a->data[1] - log_lik_a->data[0];
  psort_data[1] = log_lik_a->data[5001] - log_lik_a->data[5000];
  psort_data[2] = log_lik_a->data[10001] - log_lik_a->data[10000];
  psort_data[3] = log_lik_a->data[15001] - log_lik_a->data[15000];
  b_abs(psort_data, psort_size, tmp_data, tmp_size);
  psort_size[0] = 1;
  psort_size[1] = 4;
  psort_data[0] = log_lik_a->data[1];
  psort_data[1] = log_lik_a->data[5001];
  psort_data[2] = log_lik_a->data[10001];
  psort_data[3] = log_lik_a->data[15001];
  b_abs(psort_data, psort_size, b_v_sk_data, b_v_sk_size);
  rdivide_helper(tmp_data, tmp_size, b_v_sk_data, rk_tmp, psort_size);
  n = psort_size[1];
  if (psort_size[1] <= 2) {
    if (psort_size[1] == 1) {
      tol_nbin = rk_tmp[0];
    } else if ((rk_tmp[0] < rk_tmp[1]) || (rtIsNaN(rk_tmp[0]) && (!rtIsNaN
                 (rk_tmp[1])))) {
      tol_nbin = rk_tmp[1];
    } else {
      tol_nbin = rk_tmp[0];
    }
  } else {
    if (!rtIsNaN(rk_tmp[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= psort_size[1])) {
        if (!rtIsNaN(rk_tmp[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      tol_nbin = rk_tmp[0];
    } else {
      tol_nbin = rk_tmp[idx - 1];
      i7 = idx + 1;
      for (k = i7; k <= n; k++) {
        d4 = rk_tmp[k - 1];
        if (tol_nbin < d4) {
          tol_nbin = d4;
        }
      }
    }
  }

  /*  update for tollerance check gaussian */
  tol_gauss = 1.0;

  /*  mixture, EM */
  /*  iterates until convergence */
  i7 = p_tot->size[0] * p_tot->size[1];
  p_tot->size[0] = Vin->size[1];
  p_tot->size[1] = 4;
  emxEnsureCapacity_real_T(p_tot, i7);
  loop_ub = Vin->size[1] << 2;
  for (i7 = 0; i7 < loop_ub; i7++) {
    p_tot->data[i7] = 0.0;
  }

  emxInit_real_T(&b1, 1);
  emxInit_real_T(&psi_a, 2);
  emxInit_real_T(&psi_b, 1);
  emxInit_real_T(&dk, 2);
  emxInit_real_T(&b_b1, 2);
  emxInit_real_T(&b_a1, 2);
  emxInit_real_T(&c_a1, 1);
  emxInit_real_T(&d_a1, 1);
  emxInit_int32_T(&c_psort, 2);
  do {
    exitg2 = 0;

    /*  initialize log likelihood */
    i7 = psi_a->size[0] * psi_a->size[1];
    psi_a->size[0] = b_varargin_2;
    loop_ub = (int)nnegbin;
    psi_a->size[1] = loop_ub;
    emxEnsureCapacity_real_T(psi_a, i7);
    loop_ub_tmp = b_varargin_2 * loop_ub;
    for (i7 = 0; i7 < loop_ub_tmp; i7++) {
      psi_a->data[i7] = 0.0;
    }

    i7 = psi_b->size[0];
    psi_b->size[0] = b_varargin_2;
    emxEnsureCapacity_real_T(psi_b, i7);
    for (i7 = 0; i7 < b_varargin_2; i7++) {
      psi_b->data[i7] = 0.0;
    }

    exitg1 = false;
    while ((!exitg1) && ((tol_nbin > 0.001) && (tol_gauss > 0.001))) {
      /*  E-step */
      for (idx = 0; idx < b_varargin_2; idx++) {
        /*   Gaussian log likelihood */
        d4 = *vB_double;
        b_sqrt(&d4);
        b_K0_double[0] = *K0_double;
        b_K0_double[1] = d4;
        psi_b_tmp = -SmRG_normlike(b_K0_double, Vin->data[idx]);
        psi_b->data[idx] = psi_b_tmp;

        /*   Neg-Bin log likelihood */
        for (n = 0; n < loop_ub; n++) {
          b_K0_double[0] = mu_sk_data[n];
          b_K0_double[1] = 1.0 / rk_data[n];
          psi_a->data[idx + psi_a->size[0] * n] = -nbinlike_mu(b_K0_double,
            Vin->data[idx]);
        }

        /*  find who's bigger between [log(P(x|a))+log(P(a)] e [log(P(x|b))+log(P(b)] */
        /*  to use the logarithm property of sum: log(sum(ai))=log_b(a0)+log_b(1+sum(b^(log_b(ai)-log_b(a0)) */
        /*  with a0>a1>...>aN */
        /*  in this way it is possible to evaluate the log posterior without */
        /*  calculating the exponential: this helps for numerical stability */
        /*  (sometimes exponential blows up) */
        /*  sorting */
        b_psort->size[0] = 1;
        b_psort->size[1] = 0;
        for (n = 0; n < loop_ub; n++) {
          d4 = meana1_data[n];
          b_log(&d4);
          i7 = b_psort->size[1];
          i8 = b_psort->size[0] * b_psort->size[1];
          b_psort->size[1] = i7 + 1;
          emxEnsureCapacity_real_T(b_psort, i8);
          b_psort->data[i7] = psi_a->data[idx + psi_a->size[0] * n] + d4;
        }

        d4 = pb;
        b_log(&d4);
        i7 = b_psort->size[1];
        i8 = b_psort->size[0] * b_psort->size[1];
        b_psort->size[1] = i7 + 1;
        emxEnsureCapacity_real_T(b_psort, i8);
        b_psort->data[i7] = psi_b_tmp + d4;
        sort(b_psort, c_psort);

        /*  exponential */
        tol_nbin = 0.0;
        i7 = (int)(nnegbin + -1.0);
        for (n = 0; n < i7; n++) {
          tol_nbin += exp(b_psort->data[n + 1] - b_psort->data[0]);
        }

        /*  den */
        d4 = 1.0 + tol_nbin;
        b_log(&d4);
        tol_gauss = b_psort->data[0] + d4;

        /*  log posterior gaussian */
        d4 = pb;
        b_log(&d4);
        b->data[idx] = (psi_b_tmp + d4) - tol_gauss;

        /*  log posterior nbin */
        for (n = 0; n < loop_ub; n++) {
          d4 = meana1_data[n];
          b_log(&d4);
          a->data[idx + a->size[0] * n] = (psi_a->data[idx + psi_a->size[0] * n]
            + d4) - tol_gauss;
        }
      }

      i7 = b1->size[0];
      b1->size[0] = b->size[0];
      emxEnsureCapacity_real_T(b1, i7);
      k = b->size[0];
      for (i7 = 0; i7 < k; i7++) {
        b1->data[i7] = b->data[i7];
      }

      b_exp(b1);
      i7 = a1->size[0] * a1->size[1];
      a1->size[0] = a->size[0];
      a1->size[1] = a->size[1];
      emxEnsureCapacity_real_T(a1, i7);
      k = a->size[0] * a->size[1];
      for (i7 = 0; i7 < k; i7++) {
        a1->data[i7] = a->data[i7];
      }

      c_exp(a1);

      /*  M-step */
      /*  update Gaussian Parameters */
      i7 = b_b1->size[0] * b_b1->size[1];
      b_b1->size[0] = 1;
      b_b1->size[1] = b1->size[0];
      emxEnsureCapacity_real_T(b_b1, i7);
      k = b1->size[0];
      for (i7 = 0; i7 < k; i7++) {
        b_b1->data[i7] = b1->data[i7] * Vin->data[i7];
      }

      *K0_double = sum(b_b1) / b_sum(b1);
      i7 = b_b1->size[0] * b_b1->size[1];
      b_b1->size[0] = 1;
      b_b1->size[1] = Vin->size[1];
      emxEnsureCapacity_real_T(b_b1, i7);
      k = Vin->size[0] * Vin->size[1];
      for (i7 = 0; i7 < k; i7++) {
        b_b1->data[i7] = Vin->data[i7] - *K0_double;
      }

      power(b_b1, b_psort);
      i7 = b_b1->size[0] * b_b1->size[1];
      b_b1->size[0] = 1;
      b_b1->size[1] = b1->size[0];
      emxEnsureCapacity_real_T(b_b1, i7);
      k = b1->size[0];
      for (i7 = 0; i7 < k; i7++) {
        b_b1->data[i7] = b1->data[i7] * b_psort->data[i7];
      }

      *vB_double = sum(b_b1) / b_sum(b1);

      /*  update Neg-Bin Parameters */
      i7 = dk->size[0] * dk->size[1];
      dk->size[0] = b_varargin_2;
      dk->size[1] = loop_ub;
      emxEnsureCapacity_real_T(dk, i7);
      for (i7 = 0; i7 < loop_ub_tmp; i7++) {
        dk->data[i7] = 0.0;
      }

      for (n = 0; n < loop_ub; n++) {
        tol_gauss = minusDoubleScalarPsi(rk_data[n]);
        i7 = dk->size[0];
        tol_nbin = rk_data[n];
        i8 = b_b1->size[0] * b_b1->size[1];
        b_b1->size[0] = 1;
        b_b1->size[1] = Vin->size[1];
        emxEnsureCapacity_real_T(b_b1, i8);
        k = Vin->size[0] * Vin->size[1];
        for (i8 = 0; i8 < k; i8++) {
          b_b1->data[i8] = tol_nbin + Vin->data[i8];
        }

        psi(b_b1, b_psort);
        tol_nbin = rk_data[n];
        rk[0] = i7;
        k = rk[0];
        for (i7 = 0; i7 < k; i7++) {
          dk->data[i7 + dk->size[0] * n] = tol_nbin * (b_psort->data[i7] -
            tol_gauss);
        }
      }

      v_sk_size[0] = 1;
      v_sk_size[1] = pk_size[1];
      k = pk_size[0] * pk_size[1];
      if (0 <= k - 1) {
        memcpy(&v_sk_data[0], &pk_data[0], (unsigned int)(k * (int)sizeof(double)));
      }

      c_log(v_sk_data, v_sk_size);
      i7 = pk_size[0] * pk_size[1];
      k = i7 - 1;
      for (i7 = 0; i7 <= k; i7++) {
        v_sk_data[i7] = (1.0 - 1.0 / (1.0 - pk_data[i7])) - 1.0 / v_sk_data[i7];
      }

      /*  pk =  (betak.*sum(tauk.*dk))./sum(tauk.*(Vin'-(1-betak).*dk)); */
      for (n = 0; n < loop_ub; n++) {
        k = a1->size[0];
        i7 = c_a1->size[0];
        c_a1->size[0] = k;
        emxEnsureCapacity_real_T(c_a1, i7);
        for (i7 = 0; i7 < k; i7++) {
          c_a1->data[i7] = a1->data[i7 + a1->size[0] * n] * dk->data[i7 +
            dk->size[0] * n];
        }

        k = a1->size[0];
        i7 = d_a1->size[0];
        d_a1->size[0] = k;
        emxEnsureCapacity_real_T(d_a1, i7);
        for (i7 = 0; i7 < k; i7++) {
          d_a1->data[i7] = a1->data[i7 + a1->size[0] * n] * (Vin->data[i7] -
            (1.0 - v_sk_data[n]) * dk->data[i7 + dk->size[0] * n]);
        }

        pk_data[n] = v_sk_data[n] * b_sum(c_a1) / b_sum(d_a1);
      }

      i7 = b_a1->size[0] * b_a1->size[1];
      b_a1->size[0] = a1->size[0];
      b_a1->size[1] = a1->size[1];
      emxEnsureCapacity_real_T(b_a1, i7);
      k = a1->size[0] * a1->size[1];
      for (i7 = 0; i7 < k; i7++) {
        b_a1->data[i7] = a1->data[i7] * dk->data[i7];
      }

      c_sum(b_a1, tmp_data, tmp_size);
      c_sum(a1, b_v_sk_data, b_v_sk_size);
      rdivide_helper(tmp_data, tmp_size, b_v_sk_data, v_sk_data, v_sk_size);
      tmp_size[0] = 1;
      tmp_size[1] = pk_size[1];
      k = pk_size[0] * pk_size[1];
      if (0 <= k - 1) {
        memcpy(&tmp_data[0], &pk_data[0], (unsigned int)(k * (int)sizeof(double)));
      }

      c_log(tmp_data, tmp_size);
      b_v_sk_size[0] = 1;
      b_v_sk_size[1] = v_sk_size[1];
      k = v_sk_size[0] * v_sk_size[1];
      for (i7 = 0; i7 < k; i7++) {
        b_v_sk_data[i7] = -v_sk_data[i7];
      }

      rdivide_helper(b_v_sk_data, b_v_sk_size, tmp_data, rk_data, rk_size);
      b_pk_data.data = &pk_data[0];
      b_pk_data.size = &pk_size[0];
      b_pk_data.allocatedSize = 4;
      b_pk_data.numDimensions = 2;
      b_pk_data.canFreeData = false;
      power(&b_pk_data, b_psort);
      b_v_sk_size[0] = 1;
      b_v_sk_size[1] = rk_size[1];
      k = rk_size[0] * rk_size[1];
      for (i7 = 0; i7 < k; i7++) {
        b_v_sk_data[i7] = rk_data[i7] * (1.0 - pk_data[i7]);
      }

      rdivide_helper(b_v_sk_data, b_v_sk_size, b_psort->data, v_sk_data,
                     v_sk_size);
      mu_sk_size[0] = 1;
      mu_sk_size[1] = pk_size[1];
      k = pk_size[1];
      for (i7 = 0; i7 < k; i7++) {
        mu_sk_data[i7] = pk_data[i7] * v_sk_data[i7];
      }

      /*  today's posterior becomes tomorrow's prior :) */
      b_mean(a1, meana1_data, meana1_size);
      pb = c_mean(b1);

      /*  USEFUL IN DEBUG MODE: plots model fitting on data hist */
      /*  figure(1), */
      /*  histogram(Vin);hold on */
      /*  plot(50000*normpdf(1:max(Vin(:)),K0_double,sqrt(vB_double))) */
      /*  for inegbin = 1:nnegbin */
      /*      plot(50000*nbinpdf_mu(1:max(Vin(:)),mu_sk(inegbin),1/rk(inegbin))) */
      /*  end */
      /*  drawnow */
      /*  hold off */
      /*  update for tollerance check negbin */
      count++;
      for (n = 0; n < loop_ub; n++) {
        k = psi_a->size[0];
        i7 = b1->size[0];
        b1->size[0] = k;
        emxEnsureCapacity_real_T(b1, i7);
        for (i7 = 0; i7 < k; i7++) {
          b1->data[i7] = psi_a->data[i7 + psi_a->size[0] * n];
        }

        log_lik_a->data[count + 5000 * n] = c_mean(b1);
      }

      k = log_lik_a->size[1];
      log_lik_a_size[0] = 1;
      log_lik_a_size[1] = k;
      for (i7 = 0; i7 < k; i7++) {
        idx = count + 5000 * i7;
        psort_data[i7] = log_lik_a->data[idx] - log_lik_a->data[idx - 1];
      }

      b_abs(psort_data, log_lik_a_size, rk_tmp, psort_size);
      n = psort_size[1];
      if (psort_size[1] <= 2) {
        if (psort_size[1] == 1) {
          tol_nbin = rk_tmp[0];
        } else if ((rk_tmp[0] < rk_tmp[1]) || (rtIsNaN(rk_tmp[0]) && (!rtIsNaN
                     (rk_tmp[1])))) {
          tol_nbin = rk_tmp[1];
        } else {
          tol_nbin = rk_tmp[0];
        }
      } else {
        if (!rtIsNaN(rk_tmp[0])) {
          idx = 1;
        } else {
          idx = 0;
          k = 2;
          exitg3 = false;
          while ((!exitg3) && (k <= psort_size[1])) {
            if (!rtIsNaN(rk_tmp[k - 1])) {
              idx = k;
              exitg3 = true;
            } else {
              k++;
            }
          }
        }

        if (idx == 0) {
          tol_nbin = rk_tmp[0];
        } else {
          tol_nbin = rk_tmp[idx - 1];
          i7 = idx + 1;
          for (k = i7; k <= n; k++) {
            d4 = rk_tmp[k - 1];
            if (tol_nbin < d4) {
              tol_nbin = d4;
            }
          }
        }
      }

      /*  update for tollerance check gaussian */
      log_lik_b[count] = c_mean(psi_b);
      tol_gauss = fabs(log_lik_b[count] - log_lik_b[count - 1]) / fabs
        (log_lik_b[count]);

      /*  prevent infinite loop */
      if (count + 1 > 4999) {
        exitg1 = true;
      }
    }

    i7 = p_tot->size[0] * p_tot->size[1];
    p_tot->size[0] = a1->size[0];
    p_tot->size[1] = a1->size[1];
    emxEnsureCapacity_real_T(p_tot, i7);
    k = a1->size[0] * a1->size[1];
    for (i7 = 0; i7 < k; i7++) {
      p_tot->data[i7] = a1->data[i7];
    }

    printf("convergence reached after %d steps \n", (short)(count + 1));
    fflush(stdout);

    /*  end of EM algorithm */
    /*  now, check if some clusters are the same */
    /*  elbow method kmeans params */
    b_mu_sk_size[0] = mu_sk_size[1];
    b_mu_sk_size[1] = 2;
    k = mu_sk_size[1];
    if (0 <= k - 1) {
      memcpy(&b_mu_sk_data[0], &mu_sk_data[0], (unsigned int)(k * (int)sizeof
              (double)));
    }

    k = v_sk_size[1];
    for (i7 = 0; i7 < k; i7++) {
      b_mu_sk_data[i7 + b_mu_sk_size[0]] = v_sk_data[i7];
    }

    SmRG_kmeans_opt(b_mu_sk_data, b_mu_sk_size, nnegbin, rk_tmp, v_opt_size,
                    C_data, psort_size, psort_data, SUMD_size, &tol_gauss);

    /*  get cluster idx to merge */
    r0.f1.size[0] = loop_ub;
    if (0 <= loop_ub - 1) {
      memset(&r0.f1.data[0], 0, (unsigned int)(loop_ub * (int)sizeof(boolean_T)));
    }

    rv0[0] = r0;
    repmat(rv0, tol_gauss, r1.data, r1.size);
    structConstructorHelper(r1.data, r1.size, b_idx.data, b_idx.size);
    if (tol_gauss < nnegbin) {
      idx = -1;
      while (idx + 1 < tol_gauss) {
        idx++;
        b_idx.data[idx].to_merge.size[0] = v_opt_size[0];
        loop_ub = v_opt_size[0];
        for (i7 = 0; i7 < loop_ub; i7++) {
          b_idx.data[idx].to_merge.data[i7] = (rk_tmp[i7] == idx + 1);
        }
      }

      mus_tmp_size[0] = 1;
      mus_tmp_size_tmp = (int)tol_gauss;
      mus_tmp_size[1] = mus_tmp_size_tmp;
      if (0 <= mus_tmp_size_tmp - 1) {
        memset(&mus_tmp_data[0], 0, (unsigned int)(mus_tmp_size_tmp * (int)
                sizeof(double)));
        memset(&sigmas_tmp_data[0], 0, (unsigned int)(mus_tmp_size_tmp * (int)
                sizeof(double)));
        memset(&p_tmp_data[0], 0, (unsigned int)(mus_tmp_size_tmp * (int)sizeof
                (double)));
      }

      i7 = b_idx.size[0];
      for (count = 0; count < i7; count++) {
        /*  merge */
        idx = b_idx.data[count].to_merge.size[0] - 1;
        loop_ub_tmp = 0;
        for (n = 0; n <= idx; n++) {
          if (b_idx.data[count].to_merge.data[n]) {
            loop_ub_tmp++;
          }
        }

        k = 0;
        for (n = 0; n <= idx; n++) {
          if (b_idx.data[count].to_merge.data[n]) {
            b_tmp_data[k] = (signed char)(n + 1);
            k++;
          }
        }

        b_v_sk_size[0] = 1;
        b_v_sk_size[1] = loop_ub_tmp;
        for (i8 = 0; i8 < loop_ub_tmp; i8++) {
          b_v_sk_data[i8] = mu_sk_data[b_tmp_data[i8] - 1];
        }

        mus_tmp_data[count] = d_mean(b_v_sk_data, b_v_sk_size);
        idx = b_idx.data[count].to_merge.size[0] - 1;
        loop_ub_tmp = 0;
        for (n = 0; n <= idx; n++) {
          if (b_idx.data[count].to_merge.data[n]) {
            loop_ub_tmp++;
          }
        }

        k = 0;
        for (n = 0; n <= idx; n++) {
          if (b_idx.data[count].to_merge.data[n]) {
            c_tmp_data[k] = (signed char)(n + 1);
            k++;
          }
        }

        b_v_sk_size[0] = 1;
        b_v_sk_size[1] = loop_ub_tmp;
        for (i8 = 0; i8 < loop_ub_tmp; i8++) {
          b_v_sk_data[i8] = v_sk_data[c_tmp_data[i8] - 1];
        }

        sigmas_tmp_data[count] = d_mean(b_v_sk_data, b_v_sk_size);
        idx = b_idx.data[count].to_merge.size[0] - 1;
        loop_ub_tmp = 0;
        for (n = 0; n <= idx; n++) {
          if (b_idx.data[count].to_merge.data[n]) {
            loop_ub_tmp++;
          }
        }

        k = 0;
        for (n = 0; n <= idx; n++) {
          if (b_idx.data[count].to_merge.data[n]) {
            d_tmp_data[k] = (signed char)(n + 1);
            k++;
          }
        }

        for (i8 = 0; i8 < loop_ub_tmp; i8++) {
          b_v_sk_data[i8] = meana1_data[d_tmp_data[i8] - 1];
        }

        rk[0] = (signed char)loop_ub_tmp;
        c_v_sk_data.data = &b_v_sk_data[0];
        c_v_sk_data.size = &rk[0];
        c_v_sk_data.allocatedSize = 4;
        c_v_sk_data.numDimensions = 1;
        c_v_sk_data.canFreeData = false;
        p_tmp_data[count] = b_sum(&c_v_sk_data);
      }

      mu_sk_size[0] = 1;
      mu_sk_size[1] = mus_tmp_size[1];
      loop_ub = mus_tmp_size[1];
      if (0 <= loop_ub - 1) {
        memcpy(&mu_sk_data[0], &mus_tmp_data[0], (unsigned int)(loop_ub * (int)
                sizeof(double)));
      }

      v_sk_size[1] = mus_tmp_size_tmp;
      if (0 <= mus_tmp_size_tmp - 1) {
        memcpy(&v_sk_data[0], &sigmas_tmp_data[0], (unsigned int)
               (mus_tmp_size_tmp * (int)sizeof(double)));
      }

      meana1_size[1] = mus_tmp_size_tmp;
      if (0 <= mus_tmp_size_tmp - 1) {
        memcpy(&meana1_data[0], &p_tmp_data[0], (unsigned int)(mus_tmp_size_tmp *
                (int)sizeof(double)));
      }

      /*  update params */
      b_mus_tmp_data.data = &mus_tmp_data[0];
      b_mus_tmp_data.size = &mus_tmp_size[0];
      b_mus_tmp_data.allocatedSize = 3;
      b_mus_tmp_data.numDimensions = 2;
      b_mus_tmp_data.canFreeData = false;
      power(&b_mus_tmp_data, b_psort);
      for (i7 = 0; i7 < mus_tmp_size_tmp; i7++) {
        b_sigmas_tmp_data[i7] = sigmas_tmp_data[i7] - mus_tmp_data[i7];
      }

      rdivide_helper(b_psort->data, b_psort->size, b_sigmas_tmp_data, rk_data,
                     rk_size);
      rdivide_helper(mus_tmp_data, mus_tmp_size, sigmas_tmp_data, pk_data,
                     pk_size);
    }

    /*  if merge clusters then run again EM */
    if (meana1_size[1] == nnegbin) {
      exitg2 = 1;
    } else {
      nnegbin = tol_gauss;

      /*  set again initial conditions */
      /*  log_likelihood */
      i7 = log_lik_a->size[0] * log_lik_a->size[1];
      log_lik_a->size[0] = 5000;
      log_lik_a->size[1] = (int)tol_gauss;
      emxEnsureCapacity_real_T(log_lik_a, i7);
      loop_ub = 5000 * (int)tol_gauss;
      for (i7 = 0; i7 < loop_ub; i7++) {
        log_lik_a->data[i7] = 0.0;
      }

      i7 = (int)tol_gauss;
      for (n = 0; n < i7; n++) {
        log_lik_a->data[5000 * n] = 0.0;
        log_lik_a->data[1 + 5000 * n] = 10.0;
      }

      memset(&log_lik_b[0], 0, 5000U * sizeof(double));
      log_lik_b[0] = 0.0;
      log_lik_b[1] = 10.0;

      /*  posteriors */
      i7 = a->size[0] * a->size[1];
      a->size[0] = varargin_2;
      a->size[1] = (int)tol_gauss;
      emxEnsureCapacity_real_T(a, i7);
      loop_ub_tmp = varargin_2 * (int)tol_gauss;
      for (i7 = 0; i7 < loop_ub_tmp; i7++) {
        a->data[i7] = 0.0;
      }

      i7 = a1->size[0] * a1->size[1];
      a1->size[0] = varargin_2;
      a1->size[1] = (int)tol_gauss;
      emxEnsureCapacity_real_T(a1, i7);
      for (i7 = 0; i7 < loop_ub_tmp; i7++) {
        a1->data[i7] = 0.0;
      }

      i7 = b->size[0];
      b->size[0] = varargin_2;
      emxEnsureCapacity_real_T(b, i7);
      for (i7 = 0; i7 < varargin_2; i7++) {
        b->data[i7] = 0.0;
      }

      /*  tol   = 1e-2; */
      count = 1;

      /*  reset tol_nbin and tol_gauss */
      loop_ub = log_lik_a->size[1];
      log_lik_a_size[0] = 1;
      log_lik_a_size[1] = loop_ub;
      for (i7 = 0; i7 < loop_ub; i7++) {
        psort_data[i7] = log_lik_a->data[1 + 5000 * i7] - log_lik_a->data[5000 *
          i7];
      }

      b_abs(psort_data, log_lik_a_size, tmp_data, tmp_size);
      loop_ub = log_lik_a->size[1];
      log_lik_a_size[0] = 1;
      log_lik_a_size[1] = loop_ub;
      for (i7 = 0; i7 < loop_ub; i7++) {
        psort_data[i7] = log_lik_a->data[1 + 5000 * i7];
      }

      b_abs(psort_data, log_lik_a_size, b_v_sk_data, b_v_sk_size);
      rdivide_helper(tmp_data, tmp_size, b_v_sk_data, rk_tmp, psort_size);
      n = psort_size[1];
      if (psort_size[1] <= 2) {
        if (psort_size[1] == 1) {
          tol_nbin = rk_tmp[0];
        } else if ((rk_tmp[0] < rk_tmp[1]) || (rtIsNaN(rk_tmp[0]) && (!rtIsNaN
                     (rk_tmp[1])))) {
          tol_nbin = rk_tmp[1];
        } else {
          tol_nbin = rk_tmp[0];
        }
      } else {
        if (!rtIsNaN(rk_tmp[0])) {
          idx = 1;
        } else {
          idx = 0;
          k = 2;
          exitg1 = false;
          while ((!exitg1) && (k <= psort_size[1])) {
            if (!rtIsNaN(rk_tmp[k - 1])) {
              idx = k;
              exitg1 = true;
            } else {
              k++;
            }
          }
        }

        if (idx == 0) {
          tol_nbin = rk_tmp[0];
        } else {
          tol_nbin = rk_tmp[idx - 1];
          i7 = idx + 1;
          for (k = i7; k <= n; k++) {
            d4 = rk_tmp[k - 1];
            if (tol_nbin < d4) {
              tol_nbin = d4;
            }
          }
        }
      }

      /*  update for tollerance check gaussian */
      tol_gauss = 1.0;
    }
  } while (exitg2 == 0);

  emxFree_int32_T(&c_psort);
  emxFree_real_T(&d_a1);
  emxFree_real_T(&c_a1);
  emxFree_real_T(&b_a1);
  emxFree_real_T(&b_b1);
  emxFree_real_T(&dk);
  emxFree_real_T(&b_psort);
  emxFree_real_T(&psi_b);
  emxFree_real_T(&psi_a);
  emxFree_real_T(&log_lik_a);
  emxFree_real_T(&b1);
  emxFree_real_T(&b);
  emxFree_real_T(&a);

  /*  p_tot(i_nan)=1; */
}

/* End of code generation (SmRG_mixtureModelFitting_multmix.c) */
