/*
 * unPackDecVar_types.h
 *
 * Code generation for function 'unPackDecVar'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T {
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};
#endif /* struct_emxArray_real_T */
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /* typedef_emxArray_real_T */

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  real_T nGrid;
  real_T nState;
  real_T nControl;
  emxArray_real_T *xIdx;
  emxArray_real_T *uIdx;
} struct0_T;
#endif /* typedef_struct0_T */

/* End of code generation (unPackDecVar_types.h) */
