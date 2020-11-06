/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_SmRG_mixtureModelFitting_newPost_info.c
 *
 * Code generation for function 'SmRG_mixtureModelFitting_newPost'
 *
 */

/* Include files */
#include "_coder_SmRG_mixtureModelFitting_newPost_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *emlrtMexFcnResolvedFunctionsInfo(void);

/* Function Definitions */
static const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[5] = {
    "789ce556cb4edb40149da0809020944dab7e05860d525991268da9201091201521649cf8864c980778c610d8f009f005ac91d8f201f90a56482cf88faa192783"
    "8dc132955b57a857f2e3f8ccccb9f7cce8ca28f7bd9a4308cda0617c9e1e3e0b83eb7240cc8ebe8fa1e711e57391710aff3c0fc68fa3bc8f0b215ecdbf18e116",
    "67127a7208984de169a6c3296636938dd343402e084e8ec1f1993626d0c014ea61b0ae10ad84a827a028f55eea40eba0ee51e4764490210903df0f955f3f54af"
    "ca5fd79b7fc58f30affdd88e601db3a3f535bff36db7b4646c0970855124d0b58d326f791498148689e58ad734ea74d31cdeda1e6b49cc99f0a145714f7a2e54",
    "b903a482a5c46cdf627052e342ced1483d8729eb294470b41ecdfb8931ee52820f20a4bf97527f22567fc25fdfe15e9340a0779d52af9250afe6d3ec9fb6496d"
    "56924f1fde9877f4198c9ff4c75fdddcf954567a0f9773f92cf574fc2bbd5ecc7a6f3d779f62f4f4b9d3fc26a9d9f3e659b90b47ab5f4af30e837d41cb411eb5",
    "049da43c500cce6afd7eccfcffb5ff4e25d4a379d6c44cf5148b7acff4f752eaff6effeda7d4fb2bfb37f8d10097d9c4a8161b6bc5af9273222ccb529e0dfc32"
    "42debdd8bf38fffe541fd9becdb62fdf2f771eb3d4d3f15efbf2c7183d7d1e35bf516566e74777d55d6b1f9517164ba6bb505b41efbf2fff023c523ea0",
    "" };

  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 3120U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *epFieldName[6] = { "Name", "NumberOfInputs", "NumberOfOutputs",
    "ConstantInputs", "FullPath", "TimeStamp" };

  const char_T *propFieldName[4] = { "Version", "ResolvedFunctions",
    "EntryPoints", "CoverageInfo" };

  xEntryPoints = emlrtCreateStructMatrix(1, 1, 6, epFieldName);
  xInputs = emlrtCreateLogicalMatrix(1, 6);
  emlrtSetField(xEntryPoints, 0, "Name", emlrtMxCreateString(
    "SmRG_mixtureModelFitting_newPost"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs", emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs", emlrtMxCreateDoubleScalar
                (6.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(xEntryPoints, 0, "FullPath", emlrtMxCreateString(
    "C:\\Users\\Aleja\\Documents\\GitHub\\SmRG\\SmRG\\functions\\SmRG_mixtureModelFitting_newPost.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp", emlrtMxCreateDoubleScalar
                (738101.72680555552));
  xResult = emlrtCreateStructMatrix(1, 1, 4, propFieldName);
  emlrtSetField(xResult, 0, "Version", emlrtMxCreateString(
    "9.9.0.1467703 (R2020b)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions", (mxArray *)
                emlrtMexFcnResolvedFunctionsInfo());
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_SmRG_mixtureModelFitting_newPost_info.c) */
