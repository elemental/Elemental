/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_LEVEL1_C_H
#define EL_BLAS_LEVEL1_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* B = A^H 
   ======= */
ElError ElAdjointMatrix_i( ElConstMatrix_i A, ElMatrix_i B );
ElError ElAdjointMatrix_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElAdjointMatrix_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElAdjointMatrix_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElAdjointMatrix_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElAdjointDistMatrix_i( ElConstDistMatrix_i A, ElDistMatrix_i B );
ElError ElAdjointDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElAdjointDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElAdjointDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElAdjointDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Y := alpha A X + Y 
   ================== */
ElError ElAxpyMatrix_i( ElInt alpha, ElConstMatrix_i X, ElMatrix_i Y );
ElError ElAxpyMatrix_s( float alpha, ElConstMatrix_s X, ElMatrix_s Y );
ElError ElAxpyMatrix_d( double alpha, ElConstMatrix_d X, ElMatrix_d Y );
ElError ElAxpyMatrix_c( complex_float alpha, ElConstMatrix_c X, ElMatrix_c Y );
ElError ElAxpyMatrix_z( complex_double alpha, ElConstMatrix_z X, ElMatrix_z Y );

ElError ElAxpyDistMatrix_i
( ElInt alpha, ElConstDistMatrix_i X, ElDistMatrix_i Y );
ElError ElAxpyDistMatrix_s
( float alpha, ElConstDistMatrix_s X, ElDistMatrix_s Y );
ElError ElAxpyDistMatrix_d
( double alpha, ElConstDistMatrix_d X, ElDistMatrix_d Y );
ElError ElAxpyDistMatrix_c
( complex_float alpha, ElConstDistMatrix_c X, ElDistMatrix_c Y );
ElError ElAxpyDistMatrix_z
( complex_double alpha, ElConstDistMatrix_z X, ElDistMatrix_z Y );

/* tri(Y) := tri(alpha A X + Y)
   ============================ */
ElError ElAxpyTriangleMatrix_i
( ElUpperOrLower uplo, ElInt alpha, ElConstMatrix_i X, ElMatrix_i Y );
ElError ElAxpyTriangleMatrix_s
( ElUpperOrLower uplo, float alpha, ElConstMatrix_s X, ElMatrix_s Y );
ElError ElAxpyTriangleMatrix_d
( ElUpperOrLower uplo, double alpha, ElConstMatrix_d X, ElMatrix_d Y );
ElError ElAxpyTriangleMatrix_c
( ElUpperOrLower uplo, complex_float alpha, ElConstMatrix_c X, ElMatrix_c Y );
ElError ElAxpyTriangleMatrix_z
( ElUpperOrLower uplo, complex_double alpha, ElConstMatrix_z X, ElMatrix_z Y );

ElError ElAxpyTriangleDistMatrix_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstDistMatrix_i X, ElDistMatrix_i Y );
ElError ElAxpyTriangleDistMatrix_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s X, ElDistMatrix_s Y );
ElError ElAxpyTriangleDistMatrix_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d X, ElDistMatrix_d Y );
ElError ElAxpyTriangleDistMatrix_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c X, ElDistMatrix_c Y );
ElError ElAxpyTriangleDistMatrix_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z X, ElDistMatrix_z Y );

/* A := Conj(A) 
   ============ */
ElError ElConjugateMatrix_c( ElMatrix_c A );
ElError ElConjugateMatrix_z( ElMatrix_z A );

ElError ElConjugateDistMatrix_c( ElDistMatrix_c A );
ElError ElConjugateDistMatrix_z( ElDistMatrix_z A );

/* B = A 
   ===== */
ElError ElCopyMatrix_i( ElConstMatrix_i A, ElMatrix_i B );
ElError ElCopyMatrix_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElCopyMatrix_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElCopyMatrix_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElCopyMatrix_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElCopyDistMatrix_i( ElConstDistMatrix_i A, ElDistMatrix_i B );
ElError ElCopyDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElCopyDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElCopyDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElCopyDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* DiagonalScale 
   ============= */
ElError ElDiagonalScaleMatrix_i
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_i d, ElMatrix_i X );
ElError ElDiagonalScaleMatrix_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X );
ElError ElDiagonalScaleMatrix_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X );
ElError ElDiagonalScaleMatrix_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X );
ElError ElDiagonalScaleMatrix_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X );

ElError ElDiagonalScaleDistMatrix_i
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_i d, ElDistMatrix_i X );
ElError ElDiagonalScaleDistMatrix_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X );
ElError ElDiagonalScaleDistMatrix_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X );
ElError ElDiagonalScaleDistMatrix_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X );
ElError ElDiagonalScaleDistMatrix_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X );

/* DiagonalScaleTrapezoid
   ====================== */
ElError ElDiagonalScaleTrapezoidMatrix_i
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_i d, ElMatrix_i X, ElInt offset );
ElError ElDiagonalScaleTrapezoidMatrix_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X, ElInt offset );
ElError ElDiagonalScaleTrapezoidMatrix_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X, ElInt offset );
ElError ElDiagonalScaleTrapezoidMatrix_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X, ElInt offset );
ElError ElDiagonalScaleTrapezoidMatrix_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X, ElInt offset );

ElError ElDiagonalScaleTrapezoidDistMatrix_i
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_i d, ElDistMatrix_i X, ElInt offset );
ElError ElDiagonalScaleTrapezoidDistMatrix_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X, ElInt offset );
ElError ElDiagonalScaleTrapezoidDistMatrix_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X, ElInt offset );
ElError ElDiagonalScaleTrapezoidDistMatrix_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X, ElInt offset );
ElError ElDiagonalScaleTrapezoidDistMatrix_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X, ElInt offset );

/* DiagonalSolve 
   ============= */
ElError ElDiagonalSolveMatrix_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X );
ElError ElDiagonalSolveMatrix_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X );
ElError ElDiagonalSolveMatrix_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X );
ElError ElDiagonalSolveMatrix_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X );

ElError ElDiagonalSolveDistMatrix_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X );
ElError ElDiagonalSolveDistMatrix_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X );
ElError ElDiagonalSolveDistMatrix_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X );
ElError ElDiagonalSolveDistMatrix_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X );

/* DiagonalSolveTrapezoid
   ====================== */
ElError ElDiagonalSolveTrapezoidMatrix_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X, ElInt offset );
ElError ElDiagonalSolveTrapezoidMatrix_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X, ElInt offset );
ElError ElDiagonalSolveTrapezoidMatrix_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X, ElInt offset );
ElError ElDiagonalSolveTrapezoidMatrix_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X, ElInt offset );

ElError ElDiagonalSolveTrapezoidDistMatrix_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X, ElInt offset );
ElError ElDiagonalSolveTrapezoidDistMatrix_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X, ElInt offset );
ElError ElDiagonalSolveTrapezoidDistMatrix_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X, ElInt offset );
ElError ElDiagonalSolveTrapezoidDistMatrix_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X, ElInt offset );

/* Dot 
   === */
ElError ElDotMatrix_i
( ElConstMatrix_i A, ElConstMatrix_i B, ElInt* prod );
ElError ElDotMatrix_s
( ElConstMatrix_s A, ElConstMatrix_s B, float* prod );
ElError ElDotMatrix_d
( ElConstMatrix_d A, ElConstMatrix_d B, double* prod );
ElError ElDotMatrix_c
( ElConstMatrix_c A, ElConstMatrix_c B, complex_float* prod );
ElError ElDotMatrix_z
( ElConstMatrix_z A, ElConstMatrix_z B, complex_double* prod );

ElError ElDotDistMatrix_i
( ElConstDistMatrix_i A, ElConstDistMatrix_i B, ElInt* prod );
ElError ElDotDistMatrix_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, float* prod );
ElError ElDotDistMatrix_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, double* prod );
ElError ElDotDistMatrix_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, complex_float* prod );
ElError ElDotDistMatrix_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, complex_double* prod );

/* Dotu
   ==== */
ElError ElDotuMatrix_i
( ElConstMatrix_i A, ElConstMatrix_i B, ElInt* prod );
ElError ElDotuMatrix_s
( ElConstMatrix_s A, ElConstMatrix_s B, float* prod );
ElError ElDotuMatrix_d
( ElConstMatrix_d A, ElConstMatrix_d B, double* prod );
ElError ElDotuMatrix_c
( ElConstMatrix_c A, ElConstMatrix_c B, complex_float* prod );
ElError ElDotuMatrix_z
( ElConstMatrix_z A, ElConstMatrix_z B, complex_double* prod );

ElError ElDotuDistMatrix_i
( ElConstDistMatrix_i A, ElConstDistMatrix_i B, ElInt* prod );
ElError ElDotuDistMatrix_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, float* prod );
ElError ElDotuDistMatrix_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, double* prod );
ElError ElDotuDistMatrix_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, complex_float* prod );
ElError ElDotuDistMatrix_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, complex_double* prod );

/* EntrywiseFill
   ============= */
ElError ElEntrywiseFillMatrix_i( ElMatrix_i A, ElInt (*fill)() );
ElError ElEntrywiseFillMatrix_s( ElMatrix_s A, float (*fill)() );
ElError ElEntrywiseFillMatrix_d( ElMatrix_d A, double (*fill)() );
ElError ElEntrywiseFillMatrix_c( ElMatrix_c A, complex_float (*fill)() );
ElError ElEntrywiseFillMatrix_z( ElMatrix_z A, complex_double (*fill)() );

ElError ElEntrywiseFillDistMatrix_i
( ElDistMatrix_i A, ElInt (*fill)() );
ElError ElEntrywiseFillDistMatrix_s
( ElDistMatrix_s A, float (*fill)() );
ElError ElEntrywiseFillDistMatrix_d
( ElDistMatrix_d A, double (*fill)() );
ElError ElEntrywiseFillDistMatrix_c
( ElDistMatrix_c A, complex_float (*fill)() );
ElError ElEntrywiseFillDistMatrix_z
( ElDistMatrix_z A, complex_double (*fill)() );

/* EntrywiseMap
   ============ */
ElError ElEntrywiseMapMatrix_i
( ElMatrix_i A, ElInt (*func)(ElInt) );
ElError ElEntrywiseMapMatrix_s
( ElMatrix_s A, float (*func)(float) );
ElError ElEntrywiseMapMatrix_d
( ElMatrix_d A, double (*func)(double) );
ElError ElEntrywiseMapMatrix_c
( ElMatrix_c A, complex_float (*func)(complex_float) );
ElError ElEntrywiseMapMatrix_z
( ElMatrix_z A, complex_double (*func)(complex_double) );

ElError ElEntrywiseMapDistMatrix_i
( ElDistMatrix_i A, ElInt (*func)(ElInt) );
ElError ElEntrywiseMapDistMatrix_s
( ElDistMatrix_s A, float (*func)(float) );
ElError ElEntrywiseMapDistMatrix_d
( ElDistMatrix_d A, double (*func)(double) );
ElError ElEntrywiseMapDistMatrix_c
( ElDistMatrix_c A, complex_float (*func)(complex_float) );
ElError ElEntrywiseMapDistMatrix_z
( ElDistMatrix_z A, complex_double (*func)(complex_double) );

/* Fill
   ==== */
ElError ElFillMatrix_i( ElMatrix_i A, ElInt alpha );
ElError ElFillMatrix_s( ElMatrix_s A, float alpha );
ElError ElFillMatrix_d( ElMatrix_d A, double alpha );
ElError ElFillMatrix_c( ElMatrix_c A, complex_float alpha );
ElError ElFillMatrix_z( ElMatrix_z A, complex_double alpha );

ElError ElFillDistMatrix_i( ElDistMatrix_i A, ElInt alpha );
ElError ElFillDistMatrix_s( ElDistMatrix_s A, float alpha );
ElError ElFillDistMatrix_d( ElDistMatrix_d A, double alpha );
ElError ElFillDistMatrix_c( ElDistMatrix_c A, complex_float alpha );
ElError ElFillDistMatrix_z( ElDistMatrix_z A, complex_double alpha );

/* Hadamard
   ======== */
ElError ElHadamardMatrix_i
( ElConstMatrix_i A, ElConstMatrix_i B, ElMatrix_i C );
ElError ElHadamardMatrix_s
( ElConstMatrix_s A, ElConstMatrix_s B, ElMatrix_s C );
ElError ElHadamardMatrix_d
( ElConstMatrix_d A, ElConstMatrix_d B, ElMatrix_d C );
ElError ElHadamardMatrix_c
( ElConstMatrix_c A, ElConstMatrix_c B, ElMatrix_c C );
ElError ElHadamardMatrix_z
( ElConstMatrix_z A, ElConstMatrix_z B, ElMatrix_z C );

ElError ElHadamardDistMatrix_i
( ElConstDistMatrix_i A, ElConstDistMatrix_i B, ElDistMatrix_i C );
ElError ElHadamardDistMatrix_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, ElDistMatrix_s C );
ElError ElHadamardDistMatrix_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, ElDistMatrix_d C );
ElError ElHadamardDistMatrix_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, ElDistMatrix_c C );
ElError ElHadamardDistMatrix_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, ElDistMatrix_z C );

/* HilbertSchmidt
   ============== */
/* NOTE: This is the same as Dot */
ElError ElHilbertSchmidtMatrix_i
( ElConstMatrix_i A, ElConstMatrix_i B, ElInt* prod );
ElError ElHilbertSchmidtMatrix_s
( ElConstMatrix_s A, ElConstMatrix_s B, float* prod );
ElError ElHilbertSchmidtMatrix_d
( ElConstMatrix_d A, ElConstMatrix_d B, double* prod );
ElError ElHilbertSchmidtMatrix_c
( ElConstMatrix_c A, ElConstMatrix_c B, complex_float* prod );
ElError ElHilbertSchmidtMatrix_z
( ElConstMatrix_z A, ElConstMatrix_z B, complex_double* prod );

ElError ElHilbertSchmidtDistMatrix_i
( ElConstDistMatrix_i A, ElConstDistMatrix_i B, ElInt* prod );
ElError ElHilbertSchmidtDistMatrix_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, float* prod );
ElError ElHilbertSchmidtDistMatrix_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, double* prod );
ElError ElHilbertSchmidtDistMatrix_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, complex_float* prod );
ElError ElHilbertSchmidtDistMatrix_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, complex_double* prod );

/* IndexDependentFill
   ================== */
ElError ElIndexDependentFillMatrix_i
( ElMatrix_i A, ElInt (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillMatrix_s
( ElMatrix_s A, float (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillMatrix_d
( ElMatrix_d A, double (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillMatrix_c
( ElMatrix_c A, complex_float (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillMatrix_z
( ElMatrix_z A, complex_double (*fill)(ElInt,ElInt) );

ElError ElIndexDependentFillDistMatrix_i
( ElDistMatrix_i A, ElInt (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillDistMatrix_s
( ElDistMatrix_s A, float (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillDistMatrix_d
( ElDistMatrix_d A, double (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillDistMatrix_c
( ElDistMatrix_c A, complex_float (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillDistMatrix_z
( ElDistMatrix_z A, complex_double (*fill)(ElInt,ElInt) );

/* IndexDependentMap
   ================= */
ElError ElIndexDependentMapMatrix_i
( ElMatrix_i A, ElInt (*func)(ElInt,ElInt,ElInt) );
ElError ElIndexDependentMapMatrix_s
( ElMatrix_s A, float (*func)(ElInt,ElInt,float) );
ElError ElIndexDependentMapMatrix_d
( ElMatrix_d A, double (*func)(ElInt,ElInt,double) );
ElError ElIndexDependentMapMatrix_c
( ElMatrix_c A, complex_float (*func)(ElInt,ElInt,complex_float) );
ElError ElIndexDependentMapMatrix_z
( ElMatrix_z A, complex_double (*func)(ElInt,ElInt,complex_double) );

ElError ElIndexDependentMapDistMatrix_i
( ElDistMatrix_i A, ElInt (*func)(ElInt,ElInt,ElInt) );
ElError ElIndexDependentMapDistMatrix_s
( ElDistMatrix_s A, float (*func)(ElInt,ElInt,float) );
ElError ElIndexDependentMapDistMatrix_d
( ElDistMatrix_d A, double (*func)(ElInt,ElInt,double) );
ElError ElIndexDependentMapDistMatrix_c
( ElDistMatrix_c A, complex_float (*func)(ElInt,ElInt,complex_float) );
ElError ElIndexDependentMapDistMatrix_z
( ElDistMatrix_z A, complex_double (*func)(ElInt,ElInt,complex_double) );

/* MakeHermitian 
   ============= */
ElError ElMakeHermitianMatrix_i( ElUpperOrLower uplo, ElMatrix_i A );
ElError ElMakeHermitianMatrix_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElMakeHermitianMatrix_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElMakeHermitianMatrix_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElMakeHermitianMatrix_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElMakeHermitianDistMatrix_i( ElUpperOrLower uplo, ElDistMatrix_i A );
ElError ElMakeHermitianDistMatrix_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElMakeHermitianDistMatrix_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElMakeHermitianDistMatrix_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElMakeHermitianDistMatrix_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* MakeReal
   ======== */
ElError ElMakeRealMatrix_c( ElMatrix_c A );
ElError ElMakeRealMatrix_z( ElMatrix_z A );

ElError ElMakeRealDistMatrix_c( ElDistMatrix_c A );
ElError ElMakeRealDistMatrix_z( ElDistMatrix_z A );

/* MakeSymmetric 
   ============= */
ElError ElMakeSymmetricMatrix_i( ElUpperOrLower uplo, ElMatrix_i A );
ElError ElMakeSymmetricMatrix_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElMakeSymmetricMatrix_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElMakeSymmetricMatrix_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElMakeSymmetricMatrix_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElMakeSymmetricDistMatrix_i( ElUpperOrLower uplo, ElDistMatrix_i A );
ElError ElMakeSymmetricDistMatrix_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElMakeSymmetricDistMatrix_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElMakeSymmetricDistMatrix_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElMakeSymmetricDistMatrix_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* MakeTrapezoidal
   =============== */
ElError ElMakeTrapezoidalMatrix_i
( ElUpperOrLower uplo, ElMatrix_i A, ElInt offset );
ElError ElMakeTrapezoidalMatrix_s
( ElUpperOrLower uplo, ElMatrix_s A, ElInt offset );
ElError ElMakeTrapezoidalMatrix_d
( ElUpperOrLower uplo, ElMatrix_d A, ElInt offset );
ElError ElMakeTrapezoidalMatrix_c
( ElUpperOrLower uplo, ElMatrix_c A, ElInt offset );
ElError ElMakeTrapezoidalMatrix_z
( ElUpperOrLower uplo, ElMatrix_z A, ElInt offset );

ElError ElMakeTrapezoidalDistMatrix_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElInt offset );
ElError ElMakeTrapezoidalDistMatrix_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElInt offset );
ElError ElMakeTrapezoidalDistMatrix_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElInt offset );
ElError ElMakeTrapezoidalDistMatrix_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElInt offset );
ElError ElMakeTrapezoidalDistMatrix_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElInt offset );

/* MakeTriangular
   ============== */
ElError ElMakeTriangularMatrix_i( ElUpperOrLower uplo, ElMatrix_i A );
ElError ElMakeTriangularMatrix_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElMakeTriangularMatrix_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElMakeTriangularMatrix_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElMakeTriangularMatrix_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElMakeTriangularDistMatrix_i( ElUpperOrLower uplo, ElDistMatrix_i A );
ElError ElMakeTriangularDistMatrix_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElMakeTriangularDistMatrix_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElMakeTriangularDistMatrix_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElMakeTriangularDistMatrix_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Nrm2
   ==== */
ElError ElNrm2Matrix_s( ElConstMatrix_s A, float* gamma );
ElError ElNrm2Matrix_d( ElConstMatrix_d A, double* gamma );
ElError ElNrm2Matrix_c( ElConstMatrix_c A, float* gamma );
ElError ElNrm2Matrix_z( ElConstMatrix_z A, double* gamma );

ElError ElNrm2DistMatrix_s( ElConstDistMatrix_s A, float* gamma );
ElError ElNrm2DistMatrix_d( ElConstDistMatrix_d A, double* gamma );
ElError ElNrm2DistMatrix_c( ElConstDistMatrix_c A, float* gamma );
ElError ElNrm2DistMatrix_z( ElConstDistMatrix_z A, double* gamma );

/* Scale
   ===== */
ElError ElScaleMatrix_i( ElInt alpha, ElMatrix_i A );
ElError ElScaleMatrix_s( float alpha, ElMatrix_s A );
ElError ElScaleMatrix_d( double alpha, ElMatrix_d A );
ElError ElScaleMatrix_c( complex_float alpha, ElMatrix_c A );
ElError ElScaleMatrix_z( complex_double alpha, ElMatrix_z A );

ElError ElScaleDistMatrix_i( ElInt alpha, ElDistMatrix_i A );
ElError ElScaleDistMatrix_s( float alpha, ElDistMatrix_s A );
ElError ElScaleDistMatrix_d( double alpha, ElDistMatrix_d A );
ElError ElScaleDistMatrix_c( complex_float alpha, ElDistMatrix_c A );
ElError ElScaleDistMatrix_z( complex_double alpha, ElDistMatrix_z A );

/* ScaleTrapezoid
   ============== */
ElError ElScaleTrapezoidMatrix_i
( ElInt alpha, ElUpperOrLower uplo, ElMatrix_i A, ElInt offset );
ElError ElScaleTrapezoidMatrix_s
( float alpha, ElUpperOrLower uplo, ElMatrix_s A, ElInt offset );
ElError ElScaleTrapezoidMatrix_d
( double alpha, ElUpperOrLower uplo, ElMatrix_d A, ElInt offset );
ElError ElScaleTrapezoidMatrix_c
( complex_float alpha, ElUpperOrLower uplo, ElMatrix_c A, ElInt offset );
ElError ElScaleTrapezoidMatrix_z
( complex_double alpha, ElUpperOrLower uplo, ElMatrix_z A, ElInt offset );

ElError ElScaleTrapezoidDistMatrix_i
( ElInt alpha, ElUpperOrLower uplo, ElDistMatrix_i A, ElInt offset );
ElError ElScaleTrapezoidDistMatrix_s
( float alpha, ElUpperOrLower uplo, ElDistMatrix_s A, ElInt offset );
ElError ElScaleTrapezoidDistMatrix_d
( double alpha, ElUpperOrLower uplo, ElDistMatrix_d A, ElInt offset );
ElError ElScaleTrapezoidDistMatrix_c
( complex_float alpha, ElUpperOrLower uplo, ElDistMatrix_c A, ElInt offset );
ElError ElScaleTrapezoidDistMatrix_z
( complex_double alpha, ElUpperOrLower uplo, ElDistMatrix_z A, ElInt offset );

/* SetDiagonal
   =========== */
ElError ElSetDiagonalMatrix_i
( ElMatrix_i A, ElInt alpha, ElInt offset );
ElError ElSetDiagonalMatrix_s
( ElMatrix_s A, float alpha, ElInt offset );
ElError ElSetDiagonalMatrix_d
( ElMatrix_d A, double alpha, ElInt offset );
ElError ElSetDiagonalMatrix_c
( ElMatrix_c A, complex_float alpha, ElInt offset );
ElError ElSetDiagonalMatrix_z
( ElMatrix_z A, complex_double alpha, ElInt offset );

ElError ElSetDiagonalDistMatrix_i
( ElDistMatrix_i A, ElInt alpha, ElInt offset );
ElError ElSetDiagonalDistMatrix_s
( ElDistMatrix_s A, float alpha, ElInt offset );
ElError ElSetDiagonalDistMatrix_d
( ElDistMatrix_d A, double alpha, ElInt offset );
ElError ElSetDiagonalDistMatrix_c
( ElDistMatrix_c A, complex_float alpha, ElInt offset );
ElError ElSetDiagonalDistMatrix_z
( ElDistMatrix_z A, complex_double alpha, ElInt offset );

/* Swap
   ==== */
ElError ElSwapMatrix_i( ElOrientation orientation, ElMatrix_i X, ElMatrix_i Y );
ElError ElSwapMatrix_s( ElOrientation orientation, ElMatrix_s X, ElMatrix_s Y );
ElError ElSwapMatrix_d( ElOrientation orientation, ElMatrix_d X, ElMatrix_d Y );
ElError ElSwapMatrix_c( ElOrientation orientation, ElMatrix_c X, ElMatrix_c Y );
ElError ElSwapMatrix_z( ElOrientation orientation, ElMatrix_z X, ElMatrix_z Y );

ElError ElSwapDistMatrix_i
( ElOrientation orientation, ElDistMatrix_i X, ElDistMatrix_i Y );
ElError ElSwapDistMatrix_s
( ElOrientation orientation, ElDistMatrix_s X, ElDistMatrix_s Y );
ElError ElSwapDistMatrix_d
( ElOrientation orientation, ElDistMatrix_d X, ElDistMatrix_d Y );
ElError ElSwapDistMatrix_c
( ElOrientation orientation, ElDistMatrix_c X, ElDistMatrix_c Y );
ElError ElSwapDistMatrix_z
( ElOrientation orientation, ElDistMatrix_z X, ElDistMatrix_z Y );

ElError ElRowSwapMatrix_i( ElMatrix_i A, ElInt to, ElInt from );
ElError ElRowSwapMatrix_s( ElMatrix_s A, ElInt to, ElInt from );
ElError ElRowSwapMatrix_d( ElMatrix_d A, ElInt to, ElInt from );
ElError ElRowSwapMatrix_c( ElMatrix_c A, ElInt to, ElInt from );
ElError ElRowSwapMatrix_z( ElMatrix_z A, ElInt to, ElInt from );

ElError ElRowSwapDistMatrix_i( ElDistMatrix_i A, ElInt to, ElInt from );
ElError ElRowSwapDistMatrix_s( ElDistMatrix_s A, ElInt to, ElInt from );
ElError ElRowSwapDistMatrix_d( ElDistMatrix_d A, ElInt to, ElInt from );
ElError ElRowSwapDistMatrix_c( ElDistMatrix_c A, ElInt to, ElInt from );
ElError ElRowSwapDistMatrix_z( ElDistMatrix_z A, ElInt to, ElInt from );

ElError ElColSwapMatrix_i( ElMatrix_i A, ElInt to, ElInt from );
ElError ElColSwapMatrix_s( ElMatrix_s A, ElInt to, ElInt from );
ElError ElColSwapMatrix_d( ElMatrix_d A, ElInt to, ElInt from );
ElError ElColSwapMatrix_c( ElMatrix_c A, ElInt to, ElInt from );
ElError ElColSwapMatrix_z( ElMatrix_z A, ElInt to, ElInt from );

ElError ElColSwapDistMatrix_i( ElDistMatrix_i A, ElInt to, ElInt from );
ElError ElColSwapDistMatrix_s( ElDistMatrix_s A, ElInt to, ElInt from );
ElError ElColSwapDistMatrix_d( ElDistMatrix_d A, ElInt to, ElInt from );
ElError ElColSwapDistMatrix_c( ElDistMatrix_c A, ElInt to, ElInt from );
ElError ElColSwapDistMatrix_z( ElDistMatrix_z A, ElInt to, ElInt from );

ElError ElSymmetricSwapMatrix_i
( ElUpperOrLower uplo, ElMatrix_i A, ElInt to, ElInt from );
ElError ElSymmetricSwapMatrix_s
( ElUpperOrLower uplo, ElMatrix_s A, ElInt to, ElInt from );
ElError ElSymmetricSwapMatrix_d
( ElUpperOrLower uplo, ElMatrix_d A, ElInt to, ElInt from );
ElError ElSymmetricSwapMatrix_c
( ElUpperOrLower uplo, ElMatrix_c A, ElInt to, ElInt from );
ElError ElSymmetricSwapMatrix_z
( ElUpperOrLower uplo, ElMatrix_z A, ElInt to, ElInt from );

ElError ElSymmetricSwapDistMatrix_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElInt to, ElInt from );
ElError ElSymmetricSwapDistMatrix_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElInt to, ElInt from );
ElError ElSymmetricSwapDistMatrix_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElInt to, ElInt from );
ElError ElSymmetricSwapDistMatrix_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElInt to, ElInt from );
ElError ElSymmetricSwapDistMatrix_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElInt to, ElInt from );

ElError ElHermitianSwapMatrix_i
( ElUpperOrLower uplo, ElMatrix_i A, ElInt to, ElInt from );
ElError ElHermitianSwapMatrix_s
( ElUpperOrLower uplo, ElMatrix_s A, ElInt to, ElInt from );
ElError ElHermitianSwapMatrix_d
( ElUpperOrLower uplo, ElMatrix_d A, ElInt to, ElInt from );
ElError ElHermitianSwapMatrix_c
( ElUpperOrLower uplo, ElMatrix_c A, ElInt to, ElInt from );
ElError ElHermitianSwapMatrix_z
( ElUpperOrLower uplo, ElMatrix_z A, ElInt to, ElInt from );

ElError ElHermitianSwapDistMatrix_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElInt to, ElInt from );
ElError ElHermitianSwapDistMatrix_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElInt to, ElInt from );
ElError ElHermitianSwapDistMatrix_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElInt to, ElInt from );
ElError ElHermitianSwapDistMatrix_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElInt to, ElInt from );
ElError ElHermitianSwapDistMatrix_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElInt to, ElInt from );

/* B = A^T 
   ======= */
ElError ElTransposeMatrix_i( ElConstMatrix_i A, ElMatrix_i B );
ElError ElTransposeMatrix_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElTransposeMatrix_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElTransposeMatrix_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElTransposeMatrix_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElTransposeDistMatrix_i( ElConstDistMatrix_i A, ElDistMatrix_i B );
ElError ElTransposeDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElTransposeDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElTransposeDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElTransposeDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* UpdateDiagonal
   ============== */
ElError ElUpdateDiagonalMatrix_i
( ElMatrix_i A, ElInt alpha, ElInt offset );
ElError ElUpdateDiagonalMatrix_s
( ElMatrix_s A, float alpha, ElInt offset );
ElError ElUpdateDiagonalMatrix_d
( ElMatrix_d A, double alpha, ElInt offset );
ElError ElUpdateDiagonalMatrix_c
( ElMatrix_c A, complex_float alpha, ElInt offset );
ElError ElUpdateDiagonalMatrix_z
( ElMatrix_z A, complex_double alpha, ElInt offset );

ElError ElUpdateDiagonalDistMatrix_i
( ElDistMatrix_i A, ElInt alpha, ElInt offset );
ElError ElUpdateDiagonalDistMatrix_s
( ElDistMatrix_s A, float alpha, ElInt offset );
ElError ElUpdateDiagonalDistMatrix_d
( ElDistMatrix_d A, double alpha, ElInt offset );
ElError ElUpdateDiagonalDistMatrix_c
( ElDistMatrix_c A, complex_float alpha, ElInt offset );
ElError ElUpdateDiagonalDistMatrix_z
( ElDistMatrix_z A, complex_double alpha, ElInt offset );

/* Zero
   ==== */
ElError ElZeroMatrix_i( ElMatrix_i A );
ElError ElZeroMatrix_s( ElMatrix_s A );
ElError ElZeroMatrix_d( ElMatrix_d A );
ElError ElZeroMatrix_c( ElMatrix_c A );
ElError ElZeroMatrix_z( ElMatrix_z A );

ElError ElZeroDistMatrix_i( ElDistMatrix_i A );
ElError ElZeroDistMatrix_s( ElDistMatrix_s A );
ElError ElZeroDistMatrix_d( ElDistMatrix_d A );
ElError ElZeroDistMatrix_c( ElDistMatrix_c A );
ElError ElZeroDistMatrix_z( ElDistMatrix_z A );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_BLAS_LEVEL1_C_H */
