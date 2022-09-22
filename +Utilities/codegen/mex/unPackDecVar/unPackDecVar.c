/*
 * unPackDecVar.c
 *
 * Code generation for function 'unPackDecVar'
 *
 */

/* Include files */
#include "unPackDecVar.h"
#include "indexShapeCheck.h"
#include "rt_nonfinite.h"
#include "unPackDecVar_emxutil.h"
#include "unPackDecVar_types.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    7,              /* lineNo */
    "unPackDecVar", /* fcnName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    8,              /* lineNo */
    "unPackDecVar", /* fcnName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    11,             /* lineNo */
    "unPackDecVar", /* fcnName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    12,             /* lineNo */
    "unPackDecVar", /* fcnName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI =
    {
        29,                  /* lineNo */
        "reshapeSizeChecks", /* fcnName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI =
    {
        109,               /* lineNo */
        "computeDimsData", /* fcnName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pathName */
};

static emlrtRTEInfo emlrtRTEI =
    {
        52,                  /* lineNo */
        13,                  /* colNo */
        "reshapeSizeChecks", /* fName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo b_emlrtRTEI =
    {
        57,                  /* lineNo */
        23,                  /* colNo */
        "reshapeSizeChecks", /* fName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo c_emlrtRTEI =
    {
        59,                  /* lineNo */
        23,                  /* colNo */
        "reshapeSizeChecks", /* fName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\reshapeSizeChecks.m" /* pName */
};

static emlrtDCInfo emlrtDCI = {
    7,              /* lineNo */
    7,              /* colNo */
    "unPackDecVar", /* fName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m", /* pName */
    1    /* checkKind */
};

static emlrtBCInfo emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    7,              /* lineNo */
    7,              /* colNo */
    "z",            /* aName */
    "unPackDecVar", /* fName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m", /* pName */
    0    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = {
    8,              /* lineNo */
    7,              /* colNo */
    "unPackDecVar", /* fName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m", /* pName */
    1    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,             /* iFirst */
    -1,             /* iLast */
    8,              /* lineNo */
    7,              /* colNo */
    "z",            /* aName */
    "unPackDecVar", /* fName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m", /* pName */
    0    /* checkKind */
};

static emlrtRTEInfo e_emlrtRTEI =
    {
        58,                   /* lineNo */
        23,                   /* colNo */
        "assertValidSizeArg", /* fName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\assertValidSizeArg.m" /* pName */
};

static emlrtRTEInfo f_emlrtRTEI =
    {
        64,                   /* lineNo */
        15,                   /* colNo */
        "assertValidSizeArg", /* fName */
        "D:\\Program "
        "Files\\MATLAB\\toolbox\\eml\\eml\\+coder\\+"
        "internal\\assertValidSizeArg.m" /* pName */
};

static emlrtRTEInfo g_emlrtRTEI = {
    11,             /* lineNo */
    1,              /* colNo */
    "unPackDecVar", /* fName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m" /* pName */
};

static emlrtRTEInfo h_emlrtRTEI = {
    12,             /* lineNo */
    1,              /* colNo */
    "unPackDecVar", /* fName */
    "D:\\dev\\DynamicLaptimeSim\\v2\\4WM_PLANAR_7DOF\\+Utilities\\unPackDecVar."
    "m" /* pName */
};

/* Function Definitions */
void unPackDecVar(const emlrtStack *sp, const emxArray_real_T *z,
                  const struct0_T *pack, emxArray_real_T *x, emxArray_real_T *u)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  const real_T *z_data;
  real_T b_pack;
  real_T *x_data;
  int32_T maxdimlen;
  int32_T nx;
  boolean_T out;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  z_data = z->data;
  st.site = &emlrtRSI;
  indexShapeCheck(&st, z->size[0], *(int32_T(*)[2])pack->xIdx->size);
  maxdimlen = pack->xIdx->size[0] * pack->xIdx->size[1];
  for (nx = 0; nx < maxdimlen; nx++) {
    b_pack = pack->xIdx->data[nx];
    if (b_pack != (int32_T)muDoubleScalarFloor(b_pack)) {
      emlrtIntegerCheckR2012b(b_pack, &emlrtDCI, (emlrtCTX)sp);
    }
    if (((int32_T)b_pack < 1) || ((int32_T)b_pack > z->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)b_pack, 1, z->size[0], &emlrtBCI,
                                    (emlrtCTX)sp);
    }
  }
  st.site = &b_emlrtRSI;
  indexShapeCheck(&st, z->size[0], *(int32_T(*)[2])pack->uIdx->size);
  maxdimlen = 4 * pack->uIdx->size[1];
  for (nx = 0; nx < maxdimlen; nx++) {
    b_pack = pack->uIdx->data[nx];
    if (b_pack != (int32_T)muDoubleScalarFloor(b_pack)) {
      emlrtIntegerCheckR2012b(b_pack, &b_emlrtDCI, (emlrtCTX)sp);
    }
    if (((int32_T)b_pack < 1) || ((int32_T)b_pack > z->size[0])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)b_pack, 1, z->size[0], &b_emlrtBCI,
                                    (emlrtCTX)sp);
    }
  }
  /*  make sure x and u are returned as vectors, [nState,nTime] and
   * [nControl,nTime] */
  st.site = &c_emlrtRSI;
  nx = pack->xIdx->size[0] * pack->xIdx->size[1];
  b_st.site = &f_emlrtRSI;
  c_st.site = &g_emlrtRSI;
  if ((pack->nState != muDoubleScalarFloor(pack->nState)) ||
      muDoubleScalarIsInf(pack->nState) || (pack->nState < -2.147483648E+9) ||
      (pack->nState > 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &e_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  if (pack->nState <= 0.0) {
    b_pack = 0.0;
  } else {
    b_pack = pack->nState;
  }
  if (!(b_pack <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &f_emlrtRTEI, "Coder:MATLAB:pmaxsize",
                                  "Coder:MATLAB:pmaxsize", 0);
  }
  c_st.site = &g_emlrtRSI;
  if ((pack->nGrid != muDoubleScalarFloor(pack->nGrid)) ||
      muDoubleScalarIsInf(pack->nGrid) || (pack->nGrid < -2.147483648E+9) ||
      (pack->nGrid > 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &e_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  if (pack->nGrid <= 0.0) {
    b_pack = 0.0;
  } else {
    b_pack = pack->nGrid;
  }
  if (!(b_pack <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &f_emlrtRTEI, "Coder:MATLAB:pmaxsize",
                                  "Coder:MATLAB:pmaxsize", 0);
  }
  maxdimlen = pack->xIdx->size[0];
  if (pack->xIdx->size[1] > pack->xIdx->size[0]) {
    maxdimlen = pack->xIdx->size[1];
  }
  maxdimlen = muIntScalarMax_sint32(nx, maxdimlen);
  if ((int32_T)pack->nState > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if ((int32_T)pack->nGrid > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  out = ((int32_T)pack->nState >= 0);
  if ((!out) || ((int32_T)pack->nGrid < 0)) {
    out = false;
  }
  if (!out) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "MATLAB:checkDimCommon:nonnegativeSize",
                                  "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }
  maxdimlen = (int32_T)pack->nState * (int32_T)pack->nGrid;
  if (maxdimlen != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  nx = x->size[0] * x->size[1];
  x->size[0] = (int32_T)pack->nState;
  x->size[1] = (int32_T)pack->nGrid;
  emxEnsureCapacity_real_T(sp, x, nx, &g_emlrtRTEI);
  x_data = x->data;
  for (nx = 0; nx < maxdimlen; nx++) {
    x_data[nx] = z_data[(int32_T)pack->xIdx->data[nx] - 1];
  }
  st.site = &d_emlrtRSI;
  nx = pack->uIdx->size[1] << 2;
  b_st.site = &f_emlrtRSI;
  c_st.site = &g_emlrtRSI;
  if ((pack->nControl != muDoubleScalarFloor(pack->nControl)) ||
      muDoubleScalarIsInf(pack->nControl) ||
      (pack->nControl < -2.147483648E+9) || (pack->nControl > 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &e_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  if (pack->nControl <= 0.0) {
    b_pack = 0.0;
  } else {
    b_pack = pack->nControl;
  }
  if (!(b_pack <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &f_emlrtRTEI, "Coder:MATLAB:pmaxsize",
                                  "Coder:MATLAB:pmaxsize", 0);
  }
  c_st.site = &g_emlrtRSI;
  if ((pack->nGrid != muDoubleScalarFloor(pack->nGrid)) ||
      muDoubleScalarIsInf(pack->nGrid) || (pack->nGrid < -2.147483648E+9) ||
      (pack->nGrid > 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(
        &c_st, &e_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
        "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
  }
  if (pack->nGrid <= 0.0) {
    b_pack = 0.0;
  } else {
    b_pack = pack->nGrid;
  }
  if (!(b_pack <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &f_emlrtRTEI, "Coder:MATLAB:pmaxsize",
                                  "Coder:MATLAB:pmaxsize", 0);
  }
  maxdimlen = 4;
  if (pack->uIdx->size[1] > 4) {
    maxdimlen = pack->uIdx->size[1];
  }
  maxdimlen = muIntScalarMax_sint32(nx, maxdimlen);
  if ((int32_T)pack->nControl > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  if ((int32_T)pack->nGrid > maxdimlen) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  out = ((int32_T)pack->nControl >= 0);
  if ((!out) || ((int32_T)pack->nGrid < 0)) {
    out = false;
  }
  if (!out) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "MATLAB:checkDimCommon:nonnegativeSize",
                                  "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }
  maxdimlen = (int32_T)pack->nControl * (int32_T)pack->nGrid;
  if (maxdimlen != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  nx = u->size[0] * u->size[1];
  u->size[0] = (int32_T)pack->nControl;
  u->size[1] = (int32_T)pack->nGrid;
  emxEnsureCapacity_real_T(sp, u, nx, &h_emlrtRTEI);
  x_data = u->data;
  for (nx = 0; nx < maxdimlen; nx++) {
    x_data[nx] = z_data[(int32_T)pack->uIdx->data[nx] - 1];
  }
}

/* End of code generation (unPackDecVar.c) */
