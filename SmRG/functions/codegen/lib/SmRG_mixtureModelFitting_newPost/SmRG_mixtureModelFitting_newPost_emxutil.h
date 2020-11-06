/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SmRG_mixtureModelFitting_newPost_emxutil.h
 *
 * Code generation for function 'SmRG_mixtureModelFitting_newPost_emxutil'
 *
 */

#ifndef SMRG_MIXTUREMODELFITTING_NEWPOST_EMXUTIL_H
#define SMRG_MIXTUREMODELFITTING_NEWPOST_EMXUTIL_H

/* Include files */
#include "SmRG_mixtureModelFitting_newPost_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>
#ifdef __cplusplus

extern "C" {

#endif

  /* Function Declarations */
  extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);
  extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
  extern void emxEnsureCapacity_uint32_T(emxArray_uint32_T *emxArray, int
    oldNumel);
  extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
  extern void emxFree_real_T(emxArray_real_T **pEmxArray);
  extern void emxFree_uint32_T(emxArray_uint32_T **pEmxArray);
  extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
  extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
  extern void emxInit_uint32_T(emxArray_uint32_T **pEmxArray, int numDimensions);

#ifdef __cplusplus

}
#endif
#endif

/* End of code generation (SmRG_mixtureModelFitting_newPost_emxutil.h) */
