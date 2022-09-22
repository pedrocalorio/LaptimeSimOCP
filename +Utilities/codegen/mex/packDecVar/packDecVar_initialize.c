/*
 * packDecVar_initialize.c
 *
 * Code generation for function 'packDecVar_initialize'
 *
 */

/* Include files */
#include "packDecVar_initialize.h"
#include "_coder_packDecVar_mex.h"
#include "packDecVar_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void packDecVar_initialize(void)
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

/* End of code generation (packDecVar_initialize.c) */
