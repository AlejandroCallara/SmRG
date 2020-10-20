/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * kmeans.h
 *
 * Code generation for function 'kmeans'
 *
 */

#ifndef KMEANS_H
#define KMEANS_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "SmRG_mixtureModelFitting_multmix_types.h"

/* Function Declarations */
extern void b_local_kmeans(const double X_data[], const int X_size[2], int k,
  int idxbest_data[], int idxbest_size[1], double Cbest_data[], int Cbest_size[2],
  double varargout_1_data[], int varargout_1_size[1]);
extern void kmeans(double X_data[], int X_size[2], double kin, double
                   idxbest_data[], int idxbest_size[1], double Cbest_data[], int
                   Cbest_size[2], double varargout_1_data[], int
                   varargout_1_size[1]);

#endif

/* End of code generation (kmeans.h) */
