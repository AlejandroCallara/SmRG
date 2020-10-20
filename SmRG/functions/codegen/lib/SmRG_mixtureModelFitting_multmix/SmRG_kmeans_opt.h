/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_kmeans_opt.h
 *
 * Code generation for function 'SmRG_kmeans_opt'
 *
 */

#ifndef SMRG_KMEANS_OPT_H
#define SMRG_KMEANS_OPT_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "SmRG_mixtureModelFitting_multmix_types.h"

/* Function Declarations */
extern void SmRG_kmeans_opt(const double X_data[], const int X_size[2], double
  MAX, double IDX_data[], int IDX_size[1], double C_data[], int C_size[2],
  double SUMD_data[], int SUMD_size[1], double *K);

#endif

/* End of code generation (SmRG_kmeans_opt.h) */
