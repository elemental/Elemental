/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TYPES_C_H
#define EL_TYPES_C_H

#ifdef __cplusplus
#include <cstdbool>
#include <cstdio>
#else
#include <stdbool.h>
#include <stdio.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char ElByte;

/* If these are changes, you must make sure that they have 
   existing MPI datatypes. This is only sometimes true for 'long long' */
#ifdef EL_USE_64BIT_INTS
typedef long long int ElInt;
typedef long long unsigned ElUnsigned;
#else
typedef int ElInt;
typedef unsigned ElUnsigned;
#endif

typedef struct { ElInt  beg, end; } ElRange_i;
typedef struct { float  beg, end; } ElRange_s;
typedef struct { double beg, end; } ElRange_d;
typedef ElRange_i ElIndexRange;

static const ElInt EL_END = -100;
/* NOTE: EL_END is considered a constant by C++ but not by C */
static const ElIndexRange EL_ALL = { 0, -100 };

typedef enum
{
  EL_SUCCESS,
  EL_ALLOC_ERROR,
  EL_OUT_OF_BOUNDS_ERROR,
  EL_ARG_ERROR,           /* input argument error */
  EL_LOGIC_ERROR,         /* catch-all for logic errors */
  EL_RUNTIME_ERROR,       /* catch-all for runtime errors */
  EL_ERROR=-1             /* catch-all if the cause is unspecified */
} ElError;
const char* ElErrorString( ElError error );

#ifdef __cplusplus
} // extern "C"
#endif

#include "El/core/Element.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  ElInt value;
  ElInt index;
} ElValueInt_i;
typedef struct
{
  float value;
  ElInt index;
} ElValueInt_s;
typedef struct
{
  double value;
  ElInt index;
} ElValueInt_d;
typedef struct
{
  complex_float value;
  ElInt index;
} ElValueInt_c;
typedef struct
{
  complex_double value;
  ElInt index;
} ElValueInt_z;

typedef struct
{
  ElInt i, j;
  ElInt value;
} ElEntry_i;
typedef struct
{
  ElInt i, j;
  float value;
} ElEntry_s;
typedef struct
{
  ElInt i, j;
  double value;
} ElEntry_d;
typedef struct
{
  ElInt i, j;
  complex_float value;
} ElEntry_c;
typedef struct
{
  ElInt i, j;
  complex_double value;
} ElEntry_z;

/* product = rho * exp(kappa*n)
   where rho lies in (usually on) the unit circle and kappa is real-valued. */
typedef struct
{
    float rho;
    float kappa;
    ElInt n;
} ElSafeProduct_s; 
typedef struct
{
    double rho;
    double kappa;
    ElInt n;
} ElSafeProduct_d;
typedef struct
{
    complex_float rho;
    float kappa;
    ElInt n;
} ElSafeProduct_c;
typedef struct 
{
    complex_double rho;
    double kappa;
    ElInt n;
} ElSafeProduct_z;

/* The basic eigenvalue structure of a Hermitian matrix */
typedef struct 
{ 
    ElInt numPositive, numNegative, numZero; 
} ElInertiaType;

typedef enum 
{
    EL_UNCONJUGATED,
    EL_CONJUGATED
} ElConjugation;

typedef enum 
{
    EL_MC,   /* Col of a matrix distribution */
    EL_MD,   /* Diagonal of a matrix distribution */
    EL_MR,   /* Row of a matrix distribution */
    EL_VC,   /* Col-major vector distribution */
    EL_VR,   /* Row-major vector distribution */
    EL_STAR, /* Give to every process */
    EL_CIRC  /* Give to a single process */
} ElDist;

typedef enum
{
    EL_ELEMENT,
    EL_BLOCK
} ElDistWrap;

typedef enum 
{
    EL_FORWARD,
    EL_BACKWARD
} ElForwardOrBackward;

typedef enum
{
    EL_ROW_MAJOR,
    EL_COLUMN_MAJOR
} ElGridOrderType;

typedef enum
{
    EL_LEFT,
    EL_RIGHT
} ElLeftOrRight;

typedef enum
{
    EL_UNSORTED,
    EL_DESCENDING,
    EL_ASCENDING
} ElSortType;

typedef enum
{
    EL_ONE_NORM,           /* Operator one norm */
    EL_INFINITY_NORM,      /* Operator infinity norm */
    EL_ENTRYWISE_ONE_NORM, /* One-norm of vectorized matrix */
    EL_MAX_NORM,           /* Maximum entry-wise magnitude */
    EL_NUCLEAR_NORM,       /* One-norm of the singular values */
    EL_FROBENIUS_NORM,     /* Two-norm of the singular values */
    EL_TWO_NORM            /* Infinity-norm of the singular values */
} ElNormType;

typedef enum
{
    EL_NORMAL,
    EL_TRANSPOSE,
    EL_ADJOINT
} ElOrientation;

typedef enum
{
    EL_NON_UNIT,
    EL_UNIT
} ElUnitOrNonUnit;

typedef enum
{
    EL_LOWER,
    EL_UPPER
} ElUpperOrLower;

typedef enum
{
    EL_VERTICAL,
    EL_HORIZONTAL
} ElVerticalOrHorizontal;

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_TYPES_C_H */
