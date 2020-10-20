/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_kmeans_opt.c
 *
 * Code generation for function 'SmRG_kmeans_opt'
 *
 */

/* Include files */
#include <string.h>
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "SmRG_kmeans_opt.h"
#include "kmeans.h"
#include "sum.h"

/* Function Definitions */
void SmRG_kmeans_opt(const double X_data[], const int X_size[2], double MAX,
                     double IDX_data[], int IDX_size[1], double C_data[], int
                     C_size[2], double SUMD_data[], int SUMD_size[1], double *K)
{
  int jj;
  double D_data[4];
  int idx;
  int loop_ub;
  int b_X_size[2];
  int ii;
  double b_X_data[8];
  double Var_data[3];
  double unusedU1_data[4];
  int unusedU1_size[1];
  double unusedU2_data[8];
  int unusedU2_size[2];
  double dist_data[4];
  int dist_size[1];
  emxArray_real_T b_dist_data;
  double tmp;
  boolean_T x_data[3];
  int c_X_size[2];
  boolean_T exitg1;
  int i_data[3];
  emxArray_real_T c_dist_data;
  boolean_T guard1 = false;
  double u0;
  int n;
  boolean_T wasnan_data[4];
  boolean_T hadnans;
  int idx_data[4];

  /* %% [IDX,C,SUMD,K]=kmeans_opt(X,varargin) returns the output of the k-means */
  /* %% algorithm with the optimal number of clusters, as determined by the ELBOW */
  /* %% method. this function treats NaNs as missing data, and ignores any rows of X that */
  /* %% contain NaNs. */
  /* %% */
  /* %% [IDX]=kmeans_opt(X) returns the cluster membership for each datapoint in */
  /* %% vector X. */
  /* %% */
  /* %% [IDX]=kmeans_opt(X,MAX) returns the cluster membership for each datapoint in */
  /* %% vector X. The Elbow method will be tried from 1 to MAX number of */
  /* %% clusters (default: square root of the number of samples) */
  /* %% [IDX]=kmeans_opt(X,MAX,CUTOFF) returns the cluster membership for each datapoint in */
  /* %% vector X. The Elbow method will be tried from 1 to MAX number of */
  /* %% clusters and will choose the number which explains a fraction CUTOFF of */
  /* %% the variance (default: 0.95) */
  /* %% [IDX]=kmeans_opt(X,MAX,CUTOFF,REPEATS) returns the cluster membership for each datapoint in */
  /* %% vector X. The Elbow method will be tried from 1 to MAX number of */
  /* %% clusters and will choose the number which explains a fraction CUTOFF of */
  /* %% the variance, taking the best of REPEATS runs of k-means (default: 3). */
  /* %% [IDX,C]=kmeans_opt(X,varargin) returns in addition, the location of the */
  /* %% centroids of each cluster. */
  /* %% [IDX,C,SUMD]=kmeans_opt(X,varargin) returns in addition, the sum of */
  /* %% point-to-cluster-centroid distances. */
  /* %% [IDX,C,SUMD,K]=kmeans_opt(X,varargin) returns in addition, the number of */
  /* %% clusters. */
  /* %% sebastien.delandtsheer@uni.lu */
  /* %% sebdelandtsheer@gmail.com */
  /* %% Thomas.sauter@uni.lu */
  /* getting the number of samples */
  /*  if nargin>1, ToTest=cell2mat(varargin(1)); else, ToTest=ceil(sqrt(m)); end */
  /*  if nargin>2, Cutoff=cell2mat(varargin(2)); else, Cutoff=0.95; end */
  /*  if nargin>3, Repeats=cell2mat(varargin(3)); else, Repeats=3; end */
  jj = (int)MAX;
  if (0 <= jj - 1) {
    memset(&D_data[0], 0, (unsigned int)(jj * (int)sizeof(double)));
  }

  /* initialize the results matrix */
  for (idx = 0; idx < jj; idx++) {
    /* for each sample */
    b_X_size[0] = X_size[0];
    b_X_size[1] = 2;
    ii = X_size[0] * X_size[1];
    if (0 <= ii - 1) {
      memcpy(&b_X_data[0], &X_data[0], (unsigned int)(ii * (int)sizeof(double)));
    }

    kmeans(b_X_data, b_X_size, 1.0 + (double)idx, unusedU1_data, unusedU1_size,
           unusedU2_data, unusedU2_size, dist_data, dist_size);

    /* compute the sum of intra-cluster distances */
    b_dist_data.data = &dist_data[0];
    b_dist_data.size = &dist_size[0];
    b_dist_data.allocatedSize = 4;
    b_dist_data.numDimensions = 1;
    b_dist_data.canFreeData = false;
    tmp = b_sum(&b_dist_data);

    /* best so far */
    /* repeat the algo */
    c_X_size[0] = X_size[0];
    c_X_size[1] = 2;
    if (0 <= ii - 1) {
      memcpy(&b_X_data[0], &X_data[0], (unsigned int)(ii * (int)sizeof(double)));
    }

    kmeans(b_X_data, c_X_size, 1.0 + (double)idx, unusedU1_data, unusedU1_size,
           unusedU2_data, unusedU2_size, dist_data, dist_size);
    c_dist_data.data = &dist_data[0];
    c_dist_data.size = &dist_size[0];
    c_dist_data.allocatedSize = 4;
    c_dist_data.numDimensions = 1;
    c_dist_data.canFreeData = false;
    u0 = b_sum(&c_dist_data);
    if ((u0 < tmp) || rtIsNaN(tmp)) {
      tmp = u0;
    }

    /* repeat the algo */
    c_X_size[0] = X_size[0];
    c_X_size[1] = 2;
    if (0 <= ii - 1) {
      memcpy(&b_X_data[0], &X_data[0], (unsigned int)(ii * (int)sizeof(double)));
    }

    kmeans(b_X_data, c_X_size, 1.0 + (double)idx, unusedU1_data, unusedU1_size,
           unusedU2_data, unusedU2_size, dist_data, dist_size);
    c_dist_data.data = &dist_data[0];
    c_dist_data.size = &dist_size[0];
    c_dist_data.allocatedSize = 4;
    c_dist_data.numDimensions = 1;
    c_dist_data.canFreeData = false;
    u0 = b_sum(&c_dist_data);
    if ((u0 < tmp) || rtIsNaN(tmp)) {
      tmp = u0;
    }

    D_data[idx] = tmp;

    /* collect the best so far in the results vecor */
  }

  if (1 > (int)MAX - 1) {
    loop_ub = 0;
  } else {
    loop_ub = (int)MAX - 1;
  }

  ii = (2 <= jj);
  for (idx = 0; idx < loop_ub; idx++) {
    Var_data[idx] = D_data[idx] - D_data[ii + idx];
  }

  /* calculate %variance explained */
  idx = 2;
  if (loop_ub != 1) {
    idx = 1;
  }

  if ((1 == idx) && (loop_ub != 0) && (loop_ub != 1)) {
    for (idx = 0; idx <= loop_ub - 2; idx++) {
      Var_data[idx + 1] += Var_data[idx];
    }
  }

  /*  figure, plot(D) */
  tmp = D_data[0] - D_data[(int)MAX - 1];
  for (ii = 0; ii < loop_ub; ii++) {
    x_data[ii] = (Var_data[ii] / tmp > 0.95);
  }

  if (loop_ub == 0) {
    loop_ub = 0;
  } else {
    idx = 0;
    ii = 1;
    jj = 1;
    exitg1 = false;
    while ((!exitg1) && (jj <= 1)) {
      guard1 = false;
      if (x_data[ii - 1]) {
        idx++;
        i_data[idx - 1] = ii;
        if (idx >= loop_ub) {
          exitg1 = true;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        ii++;
        if (ii > loop_ub) {
          ii = 1;
          jj = 2;
        }
      }
    }

    if (loop_ub == 1) {
      if (idx == 0) {
        loop_ub = 0;
      }
    } else if (1 > idx) {
      loop_ub = 0;
    } else {
      loop_ub = idx;
    }
  }

  for (ii = 0; ii < loop_ub; ii++) {
    Var_data[ii] = i_data[ii];
  }

  /* find the best index */
  *K = 1.0 + Var_data[0];

  /* get the optimal number of clusters */
  unusedU2_size[0] = X_size[0];
  unusedU2_size[1] = 2;
  loop_ub = X_size[0] * X_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&unusedU2_data[0], &X_data[0], (unsigned int)(loop_ub * (int)sizeof
            (double)));
  }

  n = X_size[0] - 1;
  if (0 <= X_size[0] - 1) {
    memset(&wasnan_data[0], 0, (unsigned int)(X_size[0] * (int)sizeof(boolean_T)));
  }

  hadnans = false;
  for (jj = 0; jj <= n; jj++) {
    idx = 0;
    exitg1 = false;
    while ((!exitg1) && (idx < 2)) {
      if (rtIsNaN(X_data[jj + X_size[0] * idx])) {
        hadnans = true;
        wasnan_data[jj] = true;
        exitg1 = true;
      } else {
        idx++;
      }
    }
  }

  if (hadnans) {
    idx = X_size[0] - 1;
    loop_ub = 0;
    for (jj = 0; jj <= idx; jj++) {
      if (!wasnan_data[jj]) {
        loop_ub++;
      }
    }

    ii = 0;
    for (jj = 0; jj <= idx; jj++) {
      if (!wasnan_data[jj]) {
        idx_data[ii] = jj + 1;
        ii++;
      }
    }

    unusedU2_size[0] = loop_ub;
    unusedU2_size[1] = 2;
    for (ii = 0; ii < loop_ub; ii++) {
      unusedU2_data[ii] = X_data[idx_data[ii] - 1];
    }

    for (ii = 0; ii < loop_ub; ii++) {
      unusedU2_data[ii + loop_ub] = X_data[(idx_data[ii] + X_size[0]) - 1];
    }
  }

  b_local_kmeans(unusedU2_data, unusedU2_size, 1 + (int)Var_data[0], idx_data,
                 unusedU1_size, C_data, C_size, SUMD_data, SUMD_size);
  if (hadnans) {
    idx = -1;
    IDX_size[0] = X_size[0];
    for (jj = 0; jj <= n; jj++) {
      if (wasnan_data[jj]) {
        IDX_data[jj] = rtNaN;
      } else {
        idx++;
        IDX_data[jj] = idx_data[idx];
      }
    }
  } else {
    IDX_size[0] = unusedU1_size[0];
    loop_ub = unusedU1_size[0];
    for (ii = 0; ii < loop_ub; ii++) {
      IDX_data[ii] = idx_data[ii];
    }
  }

  /* now rerun one last time with the optimal number of clusters */
}

/* End of code generation (SmRG_kmeans_opt.c) */
