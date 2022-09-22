/*
 * packDecVar.c
 *
 * Code generation for function 'packDecVar'
 *
 */

/* Include files */
#include "packDecVar.h"
#include "packDecVar_emxutil.h"
#include "packDecVar_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI =
    {
        7,            /* lineNo */
        "packDecVar", /* fcnName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI =
    {
        8,            /* lineNo */
        "packDecVar", /* fcnName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI =
    {
        10,           /* lineNo */
        "packDecVar", /* fcnName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI =
    {
        29,                  /* lineNo */
        "reshapeSizeChecks", /* fcnName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI =
    {
        109,               /* lineNo */
        "computeDimsData", /* fcnName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pathName */
};

static emlrtRTEInfo emlrtRTEI =
    {
        59,                  /* lineNo */
        23,                  /* colNo */
        "reshapeSizeChecks", /* fName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo b_emlrtRTEI =
    {
        52,                  /* lineNo */
        13,                  /* colNo */
        "reshapeSizeChecks", /* fName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pName */
};

static emlrtECInfo emlrtECI =
    {
        -1,           /* nDims */
        19,           /* lineNo */
        1,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtECInfo b_emlrtECI =
    {
        -1,           /* nDims */
        18,           /* lineNo */
        1,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtBCInfo emlrtBCI =
    {
        -1,           /* iFirst */
        -1,           /* iLast */
        13,           /* lineNo */
        15,           /* colNo */
        "indz",       /* aName */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m", /* pName */
        0                          /* checkKind */
};

static emlrtBCInfo b_emlrtBCI =
    {
        -1,           /* iFirst */
        -1,           /* iLast */
        18,           /* lineNo */
        3,            /* colNo */
        "z",          /* aName */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m", /* pName */
        0                          /* checkKind */
};

static emlrtBCInfo c_emlrtBCI =
    {
        -1,           /* iFirst */
        -1,           /* iLast */
        19,           /* lineNo */
        3,            /* colNo */
        "z",          /* aName */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m", /* pName */
        0                          /* checkKind */
};

static emlrtBCInfo d_emlrtBCI =
    {
        -1,           /* iFirst */
        -1,           /* iLast */
        14,           /* lineNo */
        13,           /* colNo */
        "indz",       /* aName */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m", /* pName */
        0                          /* checkKind */
};

static emlrtRTEInfo c_emlrtRTEI =
    {
        58,                   /* lineNo */
        23,                   /* colNo */
        "assertValidSizeArg", /* fName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\assertValidSizeArg.m" /* pName */
};

static emlrtRTEInfo d_emlrtRTEI =
    {
        10,           /* lineNo */
        16,           /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtRTEInfo e_emlrtRTEI =
    {
        14,           /* lineNo */
        1,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtRTEInfo f_emlrtRTEI =
    {
        17,           /* lineNo */
        1,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtRTEInfo g_emlrtRTEI =
    {
        13,           /* lineNo */
        8,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtRTEInfo h_emlrtRTEI =
    {
        18,           /* lineNo */
        3,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtRTEInfo i_emlrtRTEI =
    {
        19,           /* lineNo */
        3,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtRTEInfo j_emlrtRTEI =
    {
        24,           /* lineNo */
        1,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI =
    {
        25,           /* lineNo */
        1,            /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI =
    {
        1,            /* lineNo */
        21,           /* colNo */
        "packDecVar", /* fName */
        "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+"
        "Utilities\\packDecVar.m" /* pName */
};

/* Function Definitions */
void packDecVar(const emlrtStack *sp, const emxArray_real_T *x,
                const emxArray_real_T *u, emxArray_real_T *z, struct0_T *pack)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  emxArray_int32_T *r4;
  emxArray_uint32_T *r;
  emxArray_uint32_T *r2;
  emxArray_uint32_T *uIdx;
  const real_T *u_data;
  const real_T *x_data;
  real_T b;
  real_T b_varargin_1;
  real_T varargin_1;
  real_T *z_data;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  int32_T maxdimlen;
  int32_T nx;
  int32_T *r5;
  uint32_T *r1;
  uint32_T *r3;
  uint32_T *uIdx_data;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  u_data = u->data;
  x_data = x->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtCTX)sp);
  st.site = &emlrtRSI;
  varargin_1 = (real_T)x->size[0] * (real_T)x->size[1];
  nx = x->size[0] * x->size[1];
  b_st.site = &d_emlrtRSI;
  c_st.site = &e_emlrtRSI;
  if (varargin_1 > 2.147483647E+9) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &c_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  maxdimlen = x->size[0];
  if (x->size[1] > x->size[0]) {
    maxdimlen = x->size[1];
  }
  maxdimlen = muIntScalarMax_sint32(nx, maxdimlen);
  if ((int32_T)varargin_1 > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (maxdimlen < 1) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if ((int32_T)varargin_1 != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  st.site = &b_emlrtRSI;
  b_varargin_1 = 4.0 * (real_T)x->size[1];
  nx = u->size[1] << 2;
  b_st.site = &d_emlrtRSI;
  c_st.site = &e_emlrtRSI;
  if (b_varargin_1 > 2.147483647E+9) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &c_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  maxdimlen = 4;
  if (u->size[1] > 4) {
    maxdimlen = u->size[1];
  }
  if ((int32_T)b_varargin_1 > muIntScalarMax_sint32(nx, maxdimlen)) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if ((int32_T)b_varargin_1 != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  b = (real_T)(u->size[1] << 2) + (real_T)(x->size[0] * x->size[1]);
  emxInit_uint32_T(sp, &r, &l_emlrtRTEI);
  r1 = r->data;
  if (b < 1.0) {
    r->size[0] = 1;
    r->size[1] = 0;
  } else {
    i = r->size[0] * r->size[1];
    r->size[0] = 1;
    i1 = (int32_T)b;
    r->size[1] = i1;
    emxEnsureCapacity_uint32_T(sp, r, i, &d_emlrtRTEI);
    r1 = r->data;
    loop_ub = i1 - 1;
    for (i = 0; i <= loop_ub; i++) {
      r1[i] = i + 1U;
    }
  }
  st.site = &c_emlrtRSI;
  nx = r->size[1];
  b_st.site = &d_emlrtRSI;
  c_st.site = &e_emlrtRSI;
  b = (real_T)x->size[0] + 4.0;
  if (b > 2.147483647E+9) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &c_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  maxdimlen = 1;
  if (r->size[1] > 1) {
    maxdimlen = r->size[1];
  }
  maxdimlen = muIntScalarMax_sint32(nx, maxdimlen);
  if (x->size[0] + 4 > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if (x->size[1] > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if ((x->size[0] + 4) * x->size[1] != r->size[1]) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  /*  index of state and control variables in the decVar vector */
  if (x->size[0] < 1) {
    loop_ub = 0;
  } else {
    if (x->size[0] > x->size[0] + 4) {
      emlrtDynamicBoundsCheckR2012b(x->size[0], 1, x->size[0] + 4, &emlrtBCI,
                                    (emlrtCTX)sp);
    }
    loop_ub = x->size[0];
  }
  emxInit_uint32_T(sp, &uIdx, &e_emlrtRTEI);
  maxdimlen = x->size[0] + 4;
  nx = x->size[1];
  i = uIdx->size[0] * uIdx->size[1];
  uIdx->size[0] = 4;
  uIdx->size[1] = x->size[1];
  emxEnsureCapacity_uint32_T(sp, uIdx, i, &e_emlrtRTEI);
  uIdx_data = uIdx->data;
  for (i = 0; i < nx; i++) {
    if (((int32_T)(x->size[0] + 1U) < 1) ||
        ((int32_T)(x->size[0] + 1U) > x->size[0] + 4)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(x->size[0] + 1U), 1,
                                    x->size[0] + 4, &d_emlrtBCI, (emlrtCTX)sp);
    }
    uIdx_data[4 * i] = r1[x->size[0] + maxdimlen * i];
    if (((int32_T)(x->size[0] + 2U) < 1) ||
        ((int32_T)(x->size[0] + 2U) > x->size[0] + 4)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(x->size[0] + 2U), 1,
                                    x->size[0] + 4, &d_emlrtBCI, (emlrtCTX)sp);
    }
    uIdx_data[4 * i + 1] = r1[(x->size[0] + maxdimlen * i) + 1];
    if (((int32_T)(x->size[0] + 3U) < 1) ||
        ((int32_T)(x->size[0] + 3U) > x->size[0] + 4)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(x->size[0] + 3U), 1,
                                    x->size[0] + 4, &d_emlrtBCI, (emlrtCTX)sp);
    }
    uIdx_data[4 * i + 2] = r1[(x->size[0] + maxdimlen * i) + 2];
    if (((int32_T)(x->size[0] + 4U) < 1) ||
        ((int32_T)(x->size[0] + 4U) > x->size[0] + 4)) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(x->size[0] + 4U), 1,
                                    x->size[0] + 4, &d_emlrtBCI, (emlrtCTX)sp);
    }
    uIdx_data[4 * i + 3] = r1[(x->size[0] + maxdimlen * i) + 3];
  }
  /*  decision variables are indexed so that the defects gradients appear as a
   * banded matrix */
  i = z->size[0];
  z->size[0] = (x->size[0] + 4) * x->size[1];
  emxEnsureCapacity_real_T(sp, z, i, &f_emlrtRTEI);
  z_data = z->data;
  nx = (x->size[0] + 4) * x->size[1];
  for (i = 0; i < nx; i++) {
    z_data[i] = 0.0;
  }
  emxInit_uint32_T(sp, &r2, &g_emlrtRTEI);
  maxdimlen = x->size[0] + 4;
  nx = x->size[1];
  i = r2->size[0] * r2->size[1];
  r2->size[0] = loop_ub;
  r2->size[1] = x->size[1];
  emxEnsureCapacity_uint32_T(sp, r2, i, &g_emlrtRTEI);
  r3 = r2->data;
  for (i = 0; i < nx; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      r3[i1 + r2->size[0] * i] = r1[i1 + maxdimlen * i];
    }
  }
  emxInit_int32_T(sp, &r4, &l_emlrtRTEI);
  nx = loop_ub * x->size[1];
  maxdimlen = (x->size[0] + 4) * x->size[1];
  i = r4->size[0];
  r4->size[0] = nx;
  emxEnsureCapacity_int32_T(sp, r4, i, &h_emlrtRTEI);
  r5 = r4->data;
  for (i = 0; i < nx; i++) {
    i1 = (int32_T)r3[i];
    if ((i1 < 1) || (i1 > maxdimlen)) {
      emlrtDynamicBoundsCheckR2012b(i1, 1, maxdimlen, &b_emlrtBCI,
                                    (emlrtCTX)sp);
    }
    r5[i] = i1 - 1;
  }
  emxFree_uint32_T(sp, &r2);
  maxdimlen = (int32_T)varargin_1;
  emlrtSubAssignSizeCheckR2012b(&r4->size[0], 1, &maxdimlen, 1, &b_emlrtECI,
                                (emlrtCTX)sp);
  nx = (int32_T)varargin_1;
  for (i = 0; i < nx; i++) {
    z_data[r5[i]] = x_data[i];
  }
  i = r4->size[0];
  r4->size[0] = x->size[1] << 2;
  emxEnsureCapacity_int32_T(sp, r4, i, &i_emlrtRTEI);
  r5 = r4->data;
  nx = x->size[1] << 2;
  for (i = 0; i < nx; i++) {
    i1 = (int32_T)uIdx_data[i];
    if ((i1 < 1) || (i1 > z->size[0])) {
      emlrtDynamicBoundsCheckR2012b(i1, 1, z->size[0], &c_emlrtBCI,
                                    (emlrtCTX)sp);
    }
    r5[i] = i1 - 1;
  }
  maxdimlen = (int32_T)b_varargin_1;
  emlrtSubAssignSizeCheckR2012b(&r4->size[0], 1, &maxdimlen, 1, &emlrtECI,
                                (emlrtCTX)sp);
  nx = (int32_T)b_varargin_1;
  for (i = 0; i < nx; i++) {
    z_data[r5[i]] = u_data[i];
  }
  emxFree_int32_T(sp, &r4);
  pack->nGrid = x->size[1];
  pack->nState = x->size[0];
  pack->nControl = 4.0;
  maxdimlen = x->size[0] + 4;
  nx = x->size[1];
  i = pack->xIdx->size[0] * pack->xIdx->size[1];
  pack->xIdx->size[0] = loop_ub;
  pack->xIdx->size[1] = x->size[1];
  emxEnsureCapacity_real_T(sp, pack->xIdx, i, &j_emlrtRTEI);
  for (i = 0; i < nx; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      pack->xIdx->data[i1 + pack->xIdx->size[0] * i] = r1[i1 + maxdimlen * i];
    }
  }
  emxFree_uint32_T(sp, &r);
  i = pack->uIdx->size[0] * pack->uIdx->size[1];
  pack->uIdx->size[0] = 4;
  pack->uIdx->size[1] = uIdx->size[1];
  emxEnsureCapacity_real_T(sp, pack->uIdx, i, &k_emlrtRTEI);
  loop_ub = 4 * uIdx->size[1];
  for (i = 0; i < loop_ub; i++) {
    pack->uIdx->data[i] = uIdx_data[i];
  }
  emxFree_uint32_T(sp, &uIdx);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtCTX)sp);
}

/* End of code generation (packDecVar.c) */
