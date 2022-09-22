/*
 * unPackDecVar_terminate.c
 *
 * Code generation for function 'unPackDecVar_terminate'
 *
 */

/* Include files */
#include "unPackDecVar_terminate.h"
#include "_coder_unPackDecVar_mex.h"
#include "rt_nonfinite.h"
#include "unPackDecVar_data.h"

/* Function Definitions */
void unPackDecVar_atexit(void)
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

void unPackDecVar_terminate(void)
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

/* End of code generation (unPackDecVar_terminate.c) */
