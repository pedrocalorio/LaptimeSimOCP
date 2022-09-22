/*
 * _coder_packDecVar_mex.c
 *
 * Code generation for function '_coder_packDecVar_mex'
 *
 */

/* Include files */
#include "_coder_packDecVar_mex.h"
#include "_coder_packDecVar_api.h"
#include "packDecVar_data.h"
#include "packDecVar_initialize.h"
#include "packDecVar_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&packDecVar_atexit);
  /* Module initialization. */
  packDecVar_initialize();
  /* Dispatch the entry-point. */
  packDecVar_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  packDecVar_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, (const char_T *)"windows-1252", true);
  return emlrtRootTLSGlobal;
}

void packDecVar_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
                            const mxArray *prhs[2])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[2];
  int32_T b_nlhs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        10, "packDecVar");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 10,
                        "packDecVar");
  }
  /* Call the function. */
  packDecVar_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }
  emlrtReturnArrays(b_nlhs, &plhs[0], &outputs[0]);
}

/* End of code generation (_coder_packDecVar_mex.c) */
