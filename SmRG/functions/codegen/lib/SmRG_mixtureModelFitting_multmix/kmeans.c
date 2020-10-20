/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * kmeans.c
 *
 * Code generation for function 'kmeans'
 *
 */

/* Include files */
#include <string.h>
#include <math.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "kmeans.h"
#include "rand.h"
#include "bsearch.h"

/* Function Declarations */
static void b_distfun(double D_data[], int D_size[2], const double X_data[],
                      const int X_size[2], const double C_data[], const int
                      C_size[2], const int crows_data[], int ncrows);
static void distfun(double D_data[], int D_size[2], const double X_data[], const
                    int X_size[2], const double C_data[], const int C_size[2],
                    int crows);
static int findchanged(int changed_data[], const int idx_data[], const int
  previdx_data[], const int moved_data[], const int moved_size[1], int nmoved);
static void gcentroids(double C_data[], int C_size[2], int counts_data[], const
  double X_data[], const int X_size[2], const int idx_data[], const int
  clusters_data[], int nclusters);
static void local_kmeans(const double X_data[], const int X_size[2], int k, int
  idxbest_data[], int idxbest_size[1], double Cbest_data[], int Cbest_size[2],
  double varargout_1_data[], int varargout_1_size[1]);
static void mindim2(const double D_data[], const int D_size[2], double d_data[],
                    int d_size[1], int idx_data[], int idx_size[1]);
static void simpleRandperm(int n, int k, int idx_data[], int idx_size[1]);

/* Function Definitions */
static void b_distfun(double D_data[], int D_size[2], const double X_data[],
                      const int X_size[2], const double C_data[], const int
                      C_size[2], const int crows_data[], int ncrows)
{
  int n;
  int i;
  int i12;
  int r;
  double Ccrj;
  double a;
  int i13;
  n = X_size[0] - 1;
  for (i = 0; i < ncrows; i++) {
    i12 = crows_data[i];
    for (r = 0; r <= n; r++) {
      a = X_data[r] - C_data[i12 - 1];
      D_data[r + D_size[0] * (i12 - 1)] = a * a;
    }

    Ccrj = C_data[(crows_data[i] + C_size[0]) - 1];
    for (r = 0; r <= n; r++) {
      a = X_data[r + X_size[0]] - Ccrj;
      i13 = r + D_size[0] * (i12 - 1);
      D_data[i13] += a * a;
    }
  }
}

static void distfun(double D_data[], int D_size[2], const double X_data[], const
                    int X_size[2], const double C_data[], const int C_size[2],
                    int crows)
{
  int n;
  int r;
  double Ccrj;
  double a;
  int i11;
  n = X_size[0] - 1;
  for (r = 0; r <= n; r++) {
    a = X_data[r] - C_data[crows - 1];
    D_data[r + D_size[0] * (crows - 1)] = a * a;
  }

  Ccrj = C_data[(crows + C_size[0]) - 1];
  for (r = 0; r <= n; r++) {
    a = X_data[r + X_size[0]] - Ccrj;
    i11 = r + D_size[0] * (crows - 1);
    D_data[i11] += a * a;
  }
}

static int findchanged(int changed_data[], const int idx_data[], const int
  previdx_data[], const int moved_data[], const int moved_size[1], int nmoved)
{
  int nchanged;
  int loop_ub;
  boolean_T b_data[4];
  int i15;
  loop_ub = (signed char)moved_size[0];
  if (0 <= loop_ub - 1) {
    memset(&b_data[0], 0, (unsigned int)(loop_ub * (int)sizeof(boolean_T)));
  }

  for (loop_ub = 0; loop_ub < nmoved; loop_ub++) {
    b_data[idx_data[moved_data[loop_ub] - 1] - 1] = true;
    b_data[previdx_data[moved_data[loop_ub] - 1] - 1] = true;
  }

  nchanged = 0;
  i15 = (signed char)moved_size[0];
  for (loop_ub = 0; loop_ub < i15; loop_ub++) {
    if (b_data[loop_ub]) {
      nchanged++;
      changed_data[nchanged - 1] = loop_ub + 1;
    }
  }

  return nchanged;
}

static void gcentroids(double C_data[], int C_size[2], int counts_data[], const
  double X_data[], const int X_size[2], const int idx_data[], const int
  clusters_data[], int nclusters)
{
  int n;
  int ic;
  int clic;
  int cc;
  int i14;
  int i;
  n = X_size[0];
  for (ic = 0; ic < nclusters; ic++) {
    counts_data[clusters_data[ic] - 1] = 0;
    C_data[clusters_data[ic] - 1] = rtNaN;
    C_data[(clusters_data[ic] + C_size[0]) - 1] = rtNaN;
  }

  for (ic = 0; ic < nclusters; ic++) {
    clic = clusters_data[ic] - 1;
    cc = 0;
    C_data[clic] = 0.0;
    i14 = clic + C_size[0];
    C_data[i14] = 0.0;
    for (i = 0; i < n; i++) {
      if (idx_data[i] == clic + 1) {
        cc++;
        C_data[clic] += X_data[i];
        C_data[i14] += X_data[i + X_size[0]];
      }
    }

    counts_data[clusters_data[ic] - 1] = cc;
    C_data[clic] /= (double)cc;
    C_data[i14] /= (double)cc;
  }
}

static void local_kmeans(const double X_data[], const int X_size[2], int k, int
  idxbest_data[], int idxbest_size[1], double Cbest_data[], int Cbest_size[2],
  double varargout_1_data[], int varargout_1_size[1])
{
  int n;
  double b_index;
  int pidx;
  int iempty;
  int D_size[2];
  double D_data[16];
  double d_data[4];
  int sampleDist_size[1];
  double sampleDist_data[5];
  boolean_T DNeedsComputing;
  int c;
  boolean_T exitg1;
  int crows_data[4];
  int nonEmpties_data[4];
  int b_n;
  int empties_data[4];
  int previdx_size_idx_0;
  int previdx_data[4];
  int moved_size[1];
  int moved_data[4];
  int nchanged;
  int changed_data[4];
  int iter;
  int exitg2;
  double totsumD;
  double b_d_data[4];
  int nidx_data[4];
  int nidx_size[1];
  n = X_size[0] - 1;
  b_index = b_rand();
  Cbest_size[0] = k;
  Cbest_size[1] = 2;
  pidx = k << 1;
  if (0 <= pidx - 1) {
    memset(&Cbest_data[0], 0, (unsigned int)(pidx * (int)sizeof(double)));
  }

  iempty = (int)(1.0 + floor(b_index * (double)X_size[0]));
  Cbest_data[0] = X_data[iempty - 1];
  Cbest_data[k] = X_data[(iempty + X_size[0]) - 1];
  D_size[0] = X_size[0];
  D_size[1] = k;
  pidx = X_size[0] * k;
  if (0 <= pidx - 1) {
    memset(&D_data[0], 0, (unsigned int)(pidx * (int)sizeof(double)));
  }

  distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size, 1);
  if (0 <= D_size[0] - 1) {
    memcpy(&d_data[0], &D_data[0], (unsigned int)(D_size[0] * (int)sizeof(double)));
  }

  idxbest_size[0] = X_size[0];
  pidx = X_size[0];
  for (iempty = 0; iempty < pidx; iempty++) {
    idxbest_data[iempty] = 1;
  }

  sampleDist_size[0] = X_size[0] + 1;
  if (0 <= X_size[0]) {
    memset(&sampleDist_data[0], 0, (unsigned int)((X_size[0] + 1) * (int)sizeof
            (double)));
  }

  DNeedsComputing = false;
  c = 2;
  exitg1 = false;
  while ((!exitg1) && (c <= k)) {
    b_index = 0.0;
    sampleDist_data[0] = 0.0;
    for (iempty = 0; iempty <= n; iempty++) {
      sampleDist_data[iempty + 1] = sampleDist_data[iempty] + d_data[iempty];
      b_index += d_data[iempty];
    }

    if ((b_index == 0.0) || (rtIsInf(b_index) || rtIsNaN(b_index))) {
      simpleRandperm(n + 1, (k - c) + 1, idxbest_data, idxbest_size);
      for (iempty = c; iempty <= k; iempty++) {
        Cbest_data[iempty - 1] = X_data[idxbest_data[iempty - c] - 1];
      }

      for (iempty = c; iempty <= k; iempty++) {
        Cbest_data[(iempty + k) - 1] = X_data[(idxbest_data[iempty - c] +
          X_size[0]) - 1];
      }

      DNeedsComputing = true;
      exitg1 = true;
    } else {
      pidx = sampleDist_size[0];
      for (iempty = 0; iempty < pidx; iempty++) {
        sampleDist_data[iempty] /= b_index;
      }

      pidx = b_bsearch(sampleDist_data, sampleDist_size, b_rand());
      b_index = sampleDist_data[pidx - 1];
      if (b_index < 1.0) {
        while ((pidx <= n + 1) && (sampleDist_data[pidx] <= b_index)) {
          pidx++;
        }
      } else {
        while ((pidx >= 2) && (sampleDist_data[pidx - 2] >= b_index)) {
          pidx--;
        }
      }

      Cbest_data[c - 1] = X_data[pidx - 1];
      Cbest_data[(c + k) - 1] = X_data[(pidx + X_size[0]) - 1];
      distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size, c);
      for (iempty = 0; iempty <= n; iempty++) {
        b_index = D_data[iempty + D_size[0] * (c - 1)];
        if (b_index < d_data[iempty]) {
          d_data[iempty] = b_index;
          idxbest_data[iempty] = c;
        }
      }

      c++;
    }
  }

  if (DNeedsComputing) {
    for (c = 0; c < k; c++) {
      crows_data[c] = c + 1;
    }

    b_distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size, crows_data,
              k);
    mindim2(D_data, D_size, d_data, sampleDist_size, idxbest_data, idxbest_size);
  }

  if (0 <= k - 1) {
    memset(&crows_data[0], 0, (unsigned int)(k * (int)sizeof(int)));
  }

  for (iempty = 0; iempty <= n; iempty++) {
    crows_data[idxbest_data[iempty] - 1]++;
  }

  if (0 <= k - 1) {
    memset(&nonEmpties_data[0], 0, (unsigned int)(k * (int)sizeof(int)));
  }

  b_n = X_size[0] - 1;
  if (0 <= k - 1) {
    memset(&empties_data[0], 0, (unsigned int)(k * (int)sizeof(int)));
  }

  previdx_size_idx_0 = X_size[0];
  if (0 <= X_size[0] - 1) {
    memset(&previdx_data[0], 0, (unsigned int)(X_size[0] * (int)sizeof(int)));
  }

  moved_size[0] = X_size[0];
  if (0 <= X_size[0] - 1) {
    memset(&moved_data[0], 0, (unsigned int)(X_size[0] * (int)sizeof(int)));
  }

  for (pidx = 0; pidx < k; pidx++) {
    changed_data[pidx] = pidx + 1;
  }

  nchanged = k;
  b_index = rtInf;
  iter = 0;
  do {
    exitg2 = 0;
    iter++;
    gcentroids(Cbest_data, Cbest_size, crows_data, X_data, X_size, idxbest_data,
               changed_data, nchanged);
    b_distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size,
              changed_data, nchanged);
    c = 0;
    for (pidx = 0; pidx < nchanged; pidx++) {
      if (crows_data[changed_data[pidx] - 1] == 0) {
        c++;
        empties_data[c - 1] = changed_data[pidx];
      }
    }

    if (c > 0) {
      for (pidx = 0; pidx < c; pidx++) {
        for (iempty = 0; iempty <= b_n; iempty++) {
          D_data[iempty + D_size[0] * (empties_data[pidx] - 1)] = rtNaN;
        }
      }

      nchanged -= c;
      pidx = 0;
      iempty = 0;
      while (pidx + 1 <= nchanged) {
        if (changed_data[pidx] == empties_data[iempty]) {
          iempty++;
        } else {
          changed_data[pidx - iempty] = changed_data[pidx];
        }

        pidx++;
      }
    }

    totsumD = 0.0;
    for (iempty = 0; iempty <= b_n; iempty++) {
      totsumD += D_data[iempty + D_size[0] * (idxbest_data[iempty] - 1)];
    }

    if (b_index <= totsumD) {
      idxbest_size[0] = previdx_size_idx_0;
      if (0 <= previdx_size_idx_0 - 1) {
        memcpy(&idxbest_data[0], &previdx_data[0], (unsigned int)
               (previdx_size_idx_0 * (int)sizeof(int)));
      }

      gcentroids(Cbest_data, Cbest_size, crows_data, X_data, X_size,
                 previdx_data, changed_data, nchanged);
      exitg2 = 1;
    } else if (iter >= 100) {
      exitg2 = 1;
    } else {
      previdx_size_idx_0 = idxbest_size[0];
      if (0 <= idxbest_size[0] - 1) {
        memcpy(&previdx_data[0], &idxbest_data[0], (unsigned int)(idxbest_size[0]
                * (int)sizeof(int)));
      }

      b_index = totsumD;
      mindim2(D_data, D_size, b_d_data, sampleDist_size, nidx_data, nidx_size);
      pidx = 0;
      for (iempty = 0; iempty <= b_n; iempty++) {
        if ((nidx_data[iempty] != previdx_data[iempty]) && (D_data[iempty +
             D_size[0] * (previdx_data[iempty] - 1)] > b_d_data[iempty])) {
          pidx++;
          moved_data[pidx - 1] = iempty + 1;
          idxbest_data[iempty] = nidx_data[iempty];
        }
      }

      if (pidx == 0) {
        exitg2 = 1;
      } else {
        nchanged = findchanged(changed_data, idxbest_data, previdx_data,
          moved_data, moved_size, pidx);
      }
    }
  } while (exitg2 == 0);

  pidx = 0;
  for (iempty = 0; iempty < k; iempty++) {
    if (crows_data[iempty] > 0) {
      pidx++;
      nonEmpties_data[pidx - 1] = iempty + 1;
    }
  }

  b_distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size,
            nonEmpties_data, pidx);
  for (iempty = 0; iempty <= n; iempty++) {
    d_data[iempty] = D_data[iempty + D_size[0] * (idxbest_data[iempty] - 1)];
  }

  varargout_1_size[0] = k;
  if (0 <= k - 1) {
    memset(&varargout_1_data[0], 0, (unsigned int)(k * (int)sizeof(double)));
  }

  for (iempty = 0; iempty <= n; iempty++) {
    varargout_1_data[idxbest_data[iempty] - 1] += d_data[iempty];
  }
}

static void mindim2(const double D_data[], const int D_size[2], double d_data[],
                    int d_size[1], int idx_data[], int idx_size[1])
{
  int n;
  int k;
  int loop_ub;
  int i;
  double d3;
  n = D_size[0];
  k = D_size[1];
  d_size[0] = (signed char)D_size[0];
  loop_ub = (signed char)D_size[0];
  for (i = 0; i < loop_ub; i++) {
    d_data[i] = rtInf;
  }

  idx_size[0] = D_size[0];
  loop_ub = D_size[0];
  for (i = 0; i < loop_ub; i++) {
    idx_data[i] = 1;
  }

  for (loop_ub = 0; loop_ub < k; loop_ub++) {
    for (i = 0; i < n; i++) {
      d3 = D_data[i + D_size[0] * loop_ub];
      if (d3 < d_data[i]) {
        idx_data[i] = loop_ub + 1;
        d_data[i] = d3;
      }
    }
  }
}

static void simpleRandperm(int n, int k, int idx_data[], int idx_size[1])
{
  int t;
  int i2;
  int m;
  int numer;
  double denom;
  double pt;
  double u;
  t = 0;
  idx_size[0] = n;
  if (0 <= n - 1) {
    memset(&idx_data[0], 0, (unsigned int)(n * (int)sizeof(int)));
  }

  i2 = k - 1;
  for (m = 0; m <= i2; m++) {
    numer = k - m;
    denom = n - t;
    pt = (double)numer / denom;
    u = b_rand();
    while (u > pt) {
      t++;
      denom--;
      pt += (1.0 - pt) * ((double)numer / denom);
    }

    t++;
    denom = b_rand() * (double)(m + 1);
    denom = floor(denom);
    numer = (int)denom;
    idx_data[m] = idx_data[numer];
    idx_data[numer] = t;
  }
}

void b_local_kmeans(const double X_data[], const int X_size[2], int k, int
                    idxbest_data[], int idxbest_size[1], double Cbest_data[],
                    int Cbest_size[2], double varargout_1_data[], int
                    varargout_1_size[1])
{
  int n;
  double b_index;
  int pidx;
  int cc;
  int D_size[2];
  double D_data[16];
  double d_data[4];
  int sampleDist_size[1];
  double sampleDist_data[5];
  boolean_T DNeedsComputing;
  int c;
  boolean_T exitg1;
  int i;
  int crows_data[4];
  int nonEmpties_data[4];
  int b_n;
  int empties_data[4];
  int previdx_size_idx_0;
  int previdx_data[4];
  int moved_size[1];
  int moved_data[4];
  int nchanged;
  int changed_data[4];
  double prevtotsumD;
  int iter;
  int exitg2;
  int nempty;
  double maxd;
  int from;
  double b_d_data[4];
  int nidx_data[4];
  int nidx_size[1];
  int b_i;
  n = X_size[0] - 1;
  b_index = b_rand();
  Cbest_size[0] = k;
  Cbest_size[1] = 2;
  pidx = k << 1;
  if (0 <= pidx - 1) {
    memset(&Cbest_data[0], 0, (unsigned int)(pidx * (int)sizeof(double)));
  }

  cc = (int)(1.0 + floor(b_index * (double)X_size[0]));
  Cbest_data[0] = X_data[cc - 1];
  Cbest_data[k] = X_data[(cc + X_size[0]) - 1];
  D_size[0] = X_size[0];
  D_size[1] = k;
  pidx = X_size[0] * k;
  if (0 <= pidx - 1) {
    memset(&D_data[0], 0, (unsigned int)(pidx * (int)sizeof(double)));
  }

  distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size, 1);
  if (0 <= D_size[0] - 1) {
    memcpy(&d_data[0], &D_data[0], (unsigned int)(D_size[0] * (int)sizeof(double)));
  }

  idxbest_size[0] = X_size[0];
  pidx = X_size[0];
  for (cc = 0; cc < pidx; cc++) {
    idxbest_data[cc] = 1;
  }

  sampleDist_size[0] = X_size[0] + 1;
  if (0 <= X_size[0]) {
    memset(&sampleDist_data[0], 0, (unsigned int)((X_size[0] + 1) * (int)sizeof
            (double)));
  }

  DNeedsComputing = false;
  c = 2;
  exitg1 = false;
  while ((!exitg1) && (c <= k)) {
    b_index = 0.0;
    sampleDist_data[0] = 0.0;
    for (i = 0; i <= n; i++) {
      sampleDist_data[i + 1] = sampleDist_data[i] + d_data[i];
      b_index += d_data[i];
    }

    if ((b_index == 0.0) || (rtIsInf(b_index) || rtIsNaN(b_index))) {
      simpleRandperm(n + 1, (k - c) + 1, idxbest_data, idxbest_size);
      for (i = c; i <= k; i++) {
        Cbest_data[i - 1] = X_data[idxbest_data[i - c] - 1];
      }

      for (i = c; i <= k; i++) {
        Cbest_data[(i + k) - 1] = X_data[(idxbest_data[i - c] + X_size[0]) - 1];
      }

      DNeedsComputing = true;
      exitg1 = true;
    } else {
      pidx = sampleDist_size[0];
      for (cc = 0; cc < pidx; cc++) {
        sampleDist_data[cc] /= b_index;
      }

      pidx = b_bsearch(sampleDist_data, sampleDist_size, b_rand());
      b_index = sampleDist_data[pidx - 1];
      if (b_index < 1.0) {
        while ((pidx <= n + 1) && (sampleDist_data[pidx] <= b_index)) {
          pidx++;
        }
      } else {
        while ((pidx >= 2) && (sampleDist_data[pidx - 2] >= b_index)) {
          pidx--;
        }
      }

      Cbest_data[c - 1] = X_data[pidx - 1];
      Cbest_data[(c + k) - 1] = X_data[(pidx + X_size[0]) - 1];
      distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size, c);
      for (i = 0; i <= n; i++) {
        b_index = D_data[i + D_size[0] * (c - 1)];
        if (b_index < d_data[i]) {
          d_data[i] = b_index;
          idxbest_data[i] = c;
        }
      }

      c++;
    }
  }

  if (DNeedsComputing) {
    for (c = 0; c < k; c++) {
      crows_data[c] = c + 1;
    }

    b_distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size, crows_data,
              k);
    mindim2(D_data, D_size, d_data, sampleDist_size, idxbest_data, idxbest_size);
  }

  if (0 <= k - 1) {
    memset(&crows_data[0], 0, (unsigned int)(k * (int)sizeof(int)));
  }

  for (i = 0; i <= n; i++) {
    crows_data[idxbest_data[i] - 1]++;
  }

  if (0 <= k - 1) {
    memset(&nonEmpties_data[0], 0, (unsigned int)(k * (int)sizeof(int)));
  }

  b_n = X_size[0] - 1;
  if (0 <= k - 1) {
    memset(&empties_data[0], 0, (unsigned int)(k * (int)sizeof(int)));
  }

  previdx_size_idx_0 = X_size[0];
  if (0 <= X_size[0] - 1) {
    memset(&previdx_data[0], 0, (unsigned int)(X_size[0] * (int)sizeof(int)));
  }

  moved_size[0] = X_size[0];
  if (0 <= X_size[0] - 1) {
    memset(&moved_data[0], 0, (unsigned int)(X_size[0] * (int)sizeof(int)));
  }

  for (cc = 0; cc < k; cc++) {
    changed_data[cc] = cc + 1;
  }

  nchanged = k;
  prevtotsumD = rtInf;
  iter = 0;
  do {
    exitg2 = 0;
    iter++;
    gcentroids(Cbest_data, Cbest_size, crows_data, X_data, X_size, idxbest_data,
               changed_data, nchanged);
    b_distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size,
              changed_data, nchanged);
    nempty = -1;
    for (cc = 0; cc < nchanged; cc++) {
      if (crows_data[changed_data[cc] - 1] == 0) {
        nempty++;
        empties_data[nempty] = changed_data[cc];
      }
    }

    if (nempty + 1 > 0) {
      for (i = 0; i <= nempty; i++) {
        maxd = D_data[D_size[0] * (idxbest_data[0] - 1)];
        pidx = 0;
        for (cc = 0; cc <= b_n; cc++) {
          b_index = D_data[cc + D_size[0] * (idxbest_data[cc] - 1)];
          if (b_index > maxd) {
            maxd = b_index;
            pidx = cc;
          }
        }

        from = idxbest_data[pidx] - 1;
        if (crows_data[idxbest_data[pidx] - 1] < 2) {
          cc = 0;
          exitg1 = false;
          while ((!exitg1) && (cc <= b_n)) {
            if (crows_data[cc] > 1) {
              from = cc;
              exitg1 = true;
            } else {
              cc++;
            }
          }

          cc = 0;
          exitg1 = false;
          while ((!exitg1) && (cc <= b_n)) {
            if (idxbest_data[cc] == from + 1) {
              pidx = cc;
              exitg1 = true;
            } else {
              cc++;
            }
          }
        }

        Cbest_data[empties_data[i] - 1] = X_data[pidx];
        Cbest_data[(empties_data[i] + Cbest_size[0]) - 1] = X_data[pidx +
          X_size[0]];
        crows_data[empties_data[i] - 1] = 1;
        idxbest_data[pidx] = empties_data[i];
        distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size,
                empties_data[i]);
        pidx = X_size[0];
        crows_data[from] = 0;
        Cbest_data[from] = rtNaN;
        c = from + Cbest_size[0];
        Cbest_data[c] = rtNaN;
        cc = 0;
        Cbest_data[from] = 0.0;
        Cbest_data[c] = 0.0;
        for (b_i = 0; b_i < pidx; b_i++) {
          if (idxbest_data[b_i] == from + 1) {
            cc++;
            Cbest_data[from] += X_data[b_i];
            Cbest_data[c] += X_data[b_i + X_size[0]];
          }
        }

        crows_data[from] = cc;
        Cbest_data[from] /= (double)cc;
        Cbest_data[c] /= (double)cc;
        distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size, from + 1);
        if (nchanged < k) {
          cc = 0;
          exitg1 = false;
          while ((!exitg1) && ((cc <= nchanged - 1) && (from + 1 !=
                   changed_data[cc]))) {
            if (from + 1 > changed_data[cc]) {
              for (pidx = nchanged; pidx >= cc + 1; pidx--) {
                changed_data[pidx] = changed_data[pidx - 1];
              }

              changed_data[cc] = from + 1;
              nchanged++;
              exitg1 = true;
            } else {
              cc++;
            }
          }
        }
      }
    }

    b_index = 0.0;
    for (i = 0; i <= b_n; i++) {
      b_index += D_data[i + D_size[0] * (idxbest_data[i] - 1)];
    }

    if (prevtotsumD <= b_index) {
      idxbest_size[0] = previdx_size_idx_0;
      if (0 <= previdx_size_idx_0 - 1) {
        memcpy(&idxbest_data[0], &previdx_data[0], (unsigned int)
               (previdx_size_idx_0 * (int)sizeof(int)));
      }

      gcentroids(Cbest_data, Cbest_size, crows_data, X_data, X_size,
                 previdx_data, changed_data, nchanged);
      exitg2 = 1;
    } else if (iter >= 100) {
      exitg2 = 1;
    } else {
      previdx_size_idx_0 = idxbest_size[0];
      if (0 <= idxbest_size[0] - 1) {
        memcpy(&previdx_data[0], &idxbest_data[0], (unsigned int)(idxbest_size[0]
                * (int)sizeof(int)));
      }

      prevtotsumD = b_index;
      mindim2(D_data, D_size, b_d_data, sampleDist_size, nidx_data, nidx_size);
      pidx = 0;
      for (i = 0; i <= b_n; i++) {
        if ((nidx_data[i] != previdx_data[i]) && (D_data[i + D_size[0] *
             (previdx_data[i] - 1)] > b_d_data[i])) {
          pidx++;
          moved_data[pidx - 1] = i + 1;
          idxbest_data[i] = nidx_data[i];
        }
      }

      if (pidx == 0) {
        exitg2 = 1;
      } else {
        nchanged = findchanged(changed_data, idxbest_data, previdx_data,
          moved_data, moved_size, pidx);
      }
    }
  } while (exitg2 == 0);

  pidx = 0;
  for (i = 0; i < k; i++) {
    if (crows_data[i] > 0) {
      pidx++;
      nonEmpties_data[pidx - 1] = i + 1;
    }
  }

  b_distfun(D_data, D_size, X_data, X_size, Cbest_data, Cbest_size,
            nonEmpties_data, pidx);
  for (i = 0; i <= n; i++) {
    d_data[i] = D_data[i + D_size[0] * (idxbest_data[i] - 1)];
  }

  varargout_1_size[0] = k;
  if (0 <= k - 1) {
    memset(&varargout_1_data[0], 0, (unsigned int)(k * (int)sizeof(double)));
  }

  for (i = 0; i <= n; i++) {
    varargout_1_data[idxbest_data[i] - 1] += d_data[i];
  }
}

void kmeans(double X_data[], int X_size[2], double kin, double idxbest_data[],
            int idxbest_size[1], double Cbest_data[], int Cbest_size[2], double
            varargout_1_data[], int varargout_1_size[1])
{
  int n;
  boolean_T wasnan_data[4];
  boolean_T hadnans;
  int i;
  int j;
  boolean_T exitg1;
  int idx_data[4];
  int idx_size[1];
  int trueCount;
  int partialTrueCount;
  double b_X_data[8];
  n = X_size[0];
  if (0 <= X_size[0] - 1) {
    memset(&wasnan_data[0], 0, (unsigned int)(X_size[0] * (int)sizeof(boolean_T)));
  }

  hadnans = false;
  for (i = 0; i < n; i++) {
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j < 2)) {
      if (rtIsNaN(X_data[i + X_size[0] * j])) {
        hadnans = true;
        wasnan_data[i] = true;
        exitg1 = true;
      } else {
        j++;
      }
    }
  }

  if (hadnans) {
    j = X_size[0] - 1;
    trueCount = 0;
    for (i = 0; i <= j; i++) {
      if (!wasnan_data[i]) {
        trueCount++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i <= j; i++) {
      if (!wasnan_data[i]) {
        idx_data[partialTrueCount] = i + 1;
        partialTrueCount++;
      }
    }

    for (partialTrueCount = 0; partialTrueCount < trueCount; partialTrueCount++)
    {
      b_X_data[partialTrueCount] = X_data[idx_data[partialTrueCount] - 1];
    }

    for (partialTrueCount = 0; partialTrueCount < trueCount; partialTrueCount++)
    {
      b_X_data[partialTrueCount + trueCount] = X_data[(idx_data[partialTrueCount]
        + X_size[0]) - 1];
    }

    X_size[0] = trueCount;
    X_size[1] = 2;
    j = trueCount << 1;
    if (0 <= j - 1) {
      memcpy(&X_data[0], &b_X_data[0], (unsigned int)(j * (int)sizeof(double)));
    }
  }

  local_kmeans(X_data, X_size, (int)kin, idx_data, idx_size, Cbest_data,
               Cbest_size, varargout_1_data, varargout_1_size);
  if (hadnans) {
    j = -1;
    idxbest_size[0] = n;
    for (i = 0; i < n; i++) {
      if (wasnan_data[i]) {
        idxbest_data[i] = rtNaN;
      } else {
        j++;
        idxbest_data[i] = idx_data[j];
      }
    }
  } else {
    idxbest_size[0] = idx_size[0];
    j = idx_size[0];
    for (partialTrueCount = 0; partialTrueCount < j; partialTrueCount++) {
      idxbest_data[partialTrueCount] = idx_data[partialTrueCount];
    }
  }
}

/* End of code generation (kmeans.c) */
