/*
 * unPackDecVar.h
 *
 * Code generation for function 'unPackDecVar'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "unPackDecVar_types.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void unPackDecVar(const emlrtStack *sp, const emxArray_real_T *z,
                  const struct0_T *pack, emxArray_real_T *x,
                  emxArray_real_T *u);

/* End of code generation (unPackDecVar.h) */
