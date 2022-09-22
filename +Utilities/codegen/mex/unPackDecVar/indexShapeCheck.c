/*
 * indexShapeCheck.c
 *
 * Code generation for function 'indexShapeCheck'
 *
 */

/* Include files */
#include "indexShapeCheck.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo
    e_emlrtRSI =
        {
            38,                /* lineNo */
            "indexShapeCheck", /* fcnName */
            "D:\\Program "
            "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
            "internal\\indexShapeCheck.m" /* pathName */
};

static emlrtRTEInfo
    d_emlrtRTEI =
        {
            122,           /* lineNo */
            5,             /* colNo */
            "errOrWarnIf", /* fName */
            "D:\\Program "
            "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
            "internal\\indexShapeCheck.m" /* pName */
};

/* Function Definitions */
void indexShapeCheck(const emlrtStack *sp, int32_T matrixSize,
                     const int32_T indexSize[2])
{
  emlrtStack st;
  boolean_T nonSingletonDimFound;
  st.prev = sp;
  st.tls = sp->tls;
  if (matrixSize != 1) {
    nonSingletonDimFound = (indexSize[0] != 1);
    if (indexSize[1] != 1) {
      if (nonSingletonDimFound) {
        nonSingletonDimFound = false;
      } else {
        nonSingletonDimFound = true;
      }
    }
    if (nonSingletonDimFound &&
        (((matrixSize == 1) != (indexSize[0] == 1)) || (indexSize[1] != 1))) {
      nonSingletonDimFound = true;
    } else {
      nonSingletonDimFound = false;
    }
  } else {
    nonSingletonDimFound = false;
  }
  st.site = &e_emlrtRSI;
  if (nonSingletonDimFound) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
                                  "Coder:FE:PotentialMatrixMatrix_VM",
                                  "Coder:FE:PotentialMatrixMatrix_VM", 0);
  }
}

/* End of code generation (indexShapeCheck.c) */
