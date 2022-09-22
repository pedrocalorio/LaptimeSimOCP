/*
 * unPackDecVar_initialize.c
 *
 * Code generation for function 'unPackDecVar_initialize'
 *
 */

/* Include files */
#include "unPackDecVar_initialize.h"
#include "_coder_unPackDecVar_mex.h"
#include "rt_nonfinite.h"
#include "unPackDecVar_data.h"

/* Function Definitions */
void unPackDecVar_initialize(void)
{
  static const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (unPackDecVar_initialize.c) */
