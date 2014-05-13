/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_TYPES_CINT_H
#define EL_TYPES_CINT_H

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

/* product = rho * exp(kappa*n)
   where rho lies in (usually on) the unit circle and kappa is real-valued. */
struct ElSafeProduct_s
{
    float rho;
    float kappa;
    ElInt n;
};
struct ElSafeProduct_d
{
    double rho;
    double kappa;
    ElInt n;
};
struct ElSafeProduct_c
{
    float rhoReal, rhoImag;    
    float kappa;
    ElInt n;
};
struct ElSafeProduct_z
{
    double rhoReal, rhoImag;    
    double kappa;
    ElInt n;
};

/* The basic eigenvalue structure of a Hermitian matrix */
struct ElInertiaType
{
    ElInt numPositive, numNegative, numZero;
};

enum ElConjugation
{
    EL_UNCONJUGATED,
    EL_CONJUGATED
};

enum ElDist
{
    EL_MC,   /* Col of a matrix distribution */
    EL_MD,   /* Diagonal of a matrix distribution */
    EL_MR,   /* Row of a matrix distribution */
    EL_VC,   /* Col-major vector distribution */
    EL_VR,   /* Row-major vector distribution */
    EL_STAR, /* Give to every process */
    EL_CIRC  /* Give to a single process */
};

enum ElForwardOrBackward
{
    EL_FORWARD,
    EL_BACKWARD
};

enum ElGridOrderType
{
    EL_ROW_MAJOR,
    EL_COLUMN_MAJOR
};

enum ElLeftOrRight
{
    EL_LEFT,
    EL_RIGHT
};

enum ElSortType
{
    EL_UNSORTED,
    EL_DESCENDING,
    EL_ASCENDING
};

enum ElNormType
{
    EL_ONE_NORM,           /* Operator one norm */
    EL_INFINITY_NORM,      /* Operator infinity norm */
    EL_ENTRYWISE_ONE_NORM, /* One-norm of vectorized matrix */
    EL_MAX_NORM,           /* Maximum entry-wise magnitude */
    EL_NUCLEAR_NORM,       /* One-norm of the singular values */
    EL_FROBENIUS_NORM,     /* Two-norm of the singular values */
    EL_TWO_NORM            /* Infinity-norm of the singular values */
};

enum ElOrientation
{
    EL_NORMAL,
    EL_TRANSPOSE,
    EL_ADJOINT
};

enum ElUnitOrNonUnit
{
    EL_NON_UNIT,
    EL_UNIT
};

enum ElUpperOrLower
{
    EL_LOWER,
    EL_UPPER
};

enum ElVerticalOrHorizontal
{
    EL_VERTICAL,
    EL_HORIZONTAL
};

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_TYPES_CINT_H */
