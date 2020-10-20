/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "SmRG_mixtureModelFitting_multmix.h"
#include "main.h"
#include "SmRG_mixtureModelFitting_multmix_terminate.h"
#include "SmRG_mixtureModelFitting_multmix_emxAPI.h"
#include "SmRG_mixtureModelFitting_multmix_initialize.h"

/* Function Declarations */
static emxArray_real_T *argInit_1xUnbounded_real_T(void);
static void argInit_1xd4_real_T(double result_data[], int result_size[2]);
static double argInit_real_T(void);
static void main_SmRG_mixtureModelFitting_multmix(void);

/* Function Definitions */
static emxArray_real_T *argInit_1xUnbounded_real_T(void)
{
  emxArray_real_T *result;
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreate_real_T(1, 2);

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[idx1] = argInit_real_T();
  }

  return result;
}

static void argInit_1xd4_real_T(double result_data[], int result_size[2])
{
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result_size[0] = 1;
  result_size[1] = 2;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 2; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result_data[idx1] = argInit_real_T();
  }
}

static double argInit_real_T(void)
{
  return 0.0;
}

static void main_SmRG_mixtureModelFitting_multmix(void)
{
  emxArray_real_T *p_tot;
  emxArray_real_T *a1;
  emxArray_real_T *Vin;
  double K0_double_tmp;
  double mu_sk_data[4];
  int mu_sk_size[2];
  double rk_data[4];
  int rk_size[2];
  double vB_double;
  emxInitArray_real_T(&p_tot, 2);
  emxInitArray_real_T(&a1, 2);

  /* Initialize function 'SmRG_mixtureModelFitting_multmix' input arguments. */
  /* Initialize function input argument 'Vin'. */
  Vin = argInit_1xUnbounded_real_T();
  K0_double_tmp = argInit_real_T();

  /* Initialize function input argument 'mu_sk'. */
  argInit_1xd4_real_T(mu_sk_data, mu_sk_size);

  /* Initialize function input argument 'rk'. */
  argInit_1xd4_real_T(rk_data, rk_size);

  /* Call the entry-point 'SmRG_mixtureModelFitting_multmix'. */
  vB_double = K0_double_tmp;
  SmRG_mixtureModelFitting_multmix(Vin, &K0_double_tmp, &vB_double, mu_sk_data,
    mu_sk_size, rk_data, rk_size, p_tot, a1);
  emxDestroyArray_real_T(a1);
  emxDestroyArray_real_T(p_tot);
  emxDestroyArray_real_T(Vin);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  SmRG_mixtureModelFitting_multmix_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_SmRG_mixtureModelFitting_multmix();

  /* Terminate the application.
     You do not need to do this more than one time. */
  SmRG_mixtureModelFitting_multmix_terminate();
  return 0;
}

/* End of code generation (main.c) */
