/*
 * packDecVar_terminate.c
 *
 * Code generation for function 'packDecVar_terminate'
 *
 */

/* Include files */
#include "packDecVar_terminate.h"
#include "_coder_packDecVar_mex.h"
#include "packDecVar_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void packDecVar_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void packDecVar_terminate(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (packDecVar_terminate.c) */
