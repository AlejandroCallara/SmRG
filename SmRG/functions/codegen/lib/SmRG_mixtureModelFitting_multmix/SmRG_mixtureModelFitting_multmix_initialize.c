/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_mixtureModelFitting_multmix_initialize.c
 *
 * Code generation for function 'SmRG_mixtureModelFitting_multmix_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "SmRG_mixtureModelFitting_multmix_initialize.h"
#include "eml_rand_mt19937ar_stateful.h"

/* Function Definitions */
void SmRG_mixtureModelFitting_multmix_initialize(void)
{
  rt_InitInfAndNaN(8U);
  c_eml_rand_mt19937ar_stateful_i();
}

/* End of code generation (SmRG_mixtureModelFitting_multmix_initialize.c) */
