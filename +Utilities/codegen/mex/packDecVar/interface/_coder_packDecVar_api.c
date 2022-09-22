/*
 * _coder_packDecVar_api.c
 *
 * Code generation for function '_coder_packDecVar_api'
 *
 */

/* Include files */
#include "_coder_packDecVar_api.h"
#include "packDecVar.h"
#include "packDecVar_data.h"
#include "packDecVar_emxutil.h"
#include "packDecVar_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo m_emlrtRTEI = {
    1,                       /* lineNo */
    1,                       /* colNo */
    "_coder_packDecVar_api", /* fName */
    ""                       /* pName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static const mxArray *b_emlrt_marshallOut(const struct0_T u);

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const char_T *identifier, emxArray_real_T *y);

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *x,
                             const char_T *identifier, emxArray_real_T *y);

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  e_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_emlrt_marshallOut(const struct0_T u)
{
  static const char_T *sv[5] = {"nGrid", "nState", "nControl", "xIdx", "uIdx"};
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *m;
  const mxArray *y;
  real_T *pData;
  int32_T iv[2];
  int32_T b_i;
  int32_T c_i;
  int32_T i;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 5, (const char_T **)&sv[0]));
  b_y = NULL;
  m = emlrtCreateDoubleScalar(u.nGrid);
  emlrtAssign(&b_y, m);
  emlrtSetFieldR2017b(y, 0, (const char_T *)"nGrid", b_y, 0);
  c_y = NULL;
  m = emlrtCreateDoubleScalar(u.nState);
  emlrtAssign(&c_y, m);
  emlrtSetFieldR2017b(y, 0, (const char_T *)"nState", c_y, 1);
  d_y = NULL;
  m = emlrtCreateDoubleScalar(u.nControl);
  emlrtAssign(&d_y, m);
  emlrtSetFieldR2017b(y, 0, (const char_T *)"nControl", d_y, 2);
  e_y = NULL;
  iv[0] = u.xIdx->size[0];
  iv[1] = u.xIdx->size[1];
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < u.xIdx->size[1]; b_i++) {
    for (c_i = 0; c_i < u.xIdx->size[0]; c_i++) {
      pData[i] = u.xIdx->data[c_i + u.xIdx->size[0] * b_i];
      i++;
    }
  }
  emlrtAssign(&e_y, m);
  emlrtSetFieldR2017b(y, 0, (const char_T *)"xIdx", e_y, 3);
  f_y = NULL;
  iv[0] = 4;
  iv[1] = u.uIdx->size[1];
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < u.uIdx->size[1]; b_i++) {
    pData[i] = u.uIdx->data[4 * b_i];
    pData[i + 1] = u.uIdx->data[4 * b_i + 1];
    pData[i + 2] = u.uIdx->data[4 * b_i + 2];
    pData[i + 3] = u.uIdx->data[4 * b_i + 3];
    i += 4;
  }
  emlrtAssign(&f_y, m);
  emlrtSetFieldR2017b(y, 0, (const char_T *)"uIdx", f_y, 4);
  return y;
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(u), &thisId, y);
  emlrtDestroyArray(&u);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  f_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[2] = {-1, -1};
  int32_T iv[2];
  int32_T i;
  const boolean_T bv[2] = {true, true};
  emlrtCheckVsBuiltInR2012b((emlrtCTX)sp, msgId, src, (const char_T *)"double",
                            false, 2U, (void *)&dims[0], &bv[0], &iv[0]);
  ret->allocatedSize = iv[0] * iv[1];
  i = ret->size[0] * ret->size[1];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  emxEnsureCapacity_real_T(sp, ret, i, (emlrtRTEInfo *)NULL);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *x,
                             const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(x), &thisId, y);
  emlrtDestroyArray(&x);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  static const int32_T i = 0;
  const mxArray *m;
  const mxArray *y;
  const real_T *u_data;
  u_data = u->data;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u_data[0]);
  emlrtSetDimensions((mxArray *)m, &u->size[0], 1);
  emlrtAssign(&y, m);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[2] = {4, -1};
  int32_T iv[2];
  int32_T i;
  const boolean_T bv[2] = {false, true};
  emlrtCheckVsBuiltInR2012b((emlrtCTX)sp, msgId, src, (const char_T *)"double",
                            false, 2U, (void *)&dims[0], &bv[0], &iv[0]);
  ret->allocatedSize = iv[0] * iv[1];
  i = ret->size[0] * ret->size[1];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  emxEnsureCapacity_real_T(sp, ret, i, (emlrtRTEInfo *)NULL);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

void packDecVar_api(const mxArray *const prhs[2], int32_T nlhs,
                    const mxArray *plhs[2])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  emxArray_real_T *u;
  emxArray_real_T *x;
  emxArray_real_T *z;
  struct0_T pack;
  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &x, 2, &m_emlrtRTEI);
  emxInit_real_T(&st, &u, 2, &m_emlrtRTEI);
  emxInit_real_T(&st, &z, 1, &m_emlrtRTEI);
  emxInitStruct_struct0_T(&st, &pack, &m_emlrtRTEI);
  /* Marshall function inputs */
  x->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "x", x);
  u->canFreeData = false;
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "u", u);
  /* Invoke the target function */
  packDecVar(&st, x, u, z, &pack);
  /* Marshall function outputs */
  z->canFreeData = false;
  plhs[0] = emlrt_marshallOut(z);
  emxFree_real_T(&st, &z);
  emxFree_real_T(&st, &u);
  emxFree_real_T(&st, &x);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(pack);
  }
  emxFreeStruct_struct0_T(&st, &pack);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_packDecVar_api.c) */
