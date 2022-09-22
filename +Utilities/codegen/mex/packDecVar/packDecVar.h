/*
 * packDecVar.h
 *
 * Code generation for function 'packDecVar'
 *
 */

#pragma once

/* Include files */
#include "packDecVar_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void packDecVar(const emlrtStack *sp, const emxArray_real_T *x,
                const emxArray_real_T *u, emxArray_real_T *z, struct0_T *pack);

/* End of code generation (packDecVar.h) */
