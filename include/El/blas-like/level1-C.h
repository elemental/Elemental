/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
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
ElError ElAdjoint_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElAdjoint_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElAdjointDist_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElAdjointDist_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Y := alpha A X + Y 
   ================== */
ElError ElAxpy_i( ElInt alpha, ElConstMatrix_i X, ElMatrix_i Y );
ElError ElAxpy_s( float alpha, ElConstMatrix_s X, ElMatrix_s Y );
ElError ElAxpy_d( double alpha, ElConstMatrix_d X, ElMatrix_d Y );
ElError ElAxpy_c( complex_float alpha, ElConstMatrix_c X, ElMatrix_c Y );
ElError ElAxpy_z( complex_double alpha, ElConstMatrix_z X, ElMatrix_z Y );

ElError ElAxpyDist_i
( ElInt alpha, ElConstDistMatrix_i X, ElDistMatrix_i Y );
ElError ElAxpyDist_s
( float alpha, ElConstDistMatrix_s X, ElDistMatrix_s Y );
ElError ElAxpyDist_d
( double alpha, ElConstDistMatrix_d X, ElDistMatrix_d Y );
ElError ElAxpyDist_c
( complex_float alpha, ElConstDistMatrix_c X, ElDistMatrix_c Y );
ElError ElAxpyDist_z
( complex_double alpha, ElConstDistMatrix_z X, ElDistMatrix_z Y );

/* tri(Y) := tri(alpha A X + Y)
   ============================ */
ElError ElAxpyTriangle_i
( ElUpperOrLower uplo, ElInt alpha, ElConstMatrix_i X, ElMatrix_i Y );
ElError ElAxpyTriangle_s
( ElUpperOrLower uplo, float alpha, ElConstMatrix_s X, ElMatrix_s Y );
ElError ElAxpyTriangle_d
( ElUpperOrLower uplo, double alpha, ElConstMatrix_d X, ElMatrix_d Y );
ElError ElAxpyTriangle_c
( ElUpperOrLower uplo, complex_float alpha, ElConstMatrix_c X, ElMatrix_c Y );
ElError ElAxpyTriangle_z
( ElUpperOrLower uplo, complex_double alpha, ElConstMatrix_z X, ElMatrix_z Y );

ElError ElAxpyTriangleDist_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstDistMatrix_i X, ElDistMatrix_i Y );
ElError ElAxpyTriangleDist_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s X, ElDistMatrix_s Y );
ElError ElAxpyTriangleDist_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d X, ElDistMatrix_d Y );
ElError ElAxpyTriangleDist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c X, ElDistMatrix_c Y );
ElError ElAxpyTriangleDist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z X, ElDistMatrix_z Y );

/* A := Conj(A) 
   ============ */
ElError ElConjugate_c( ElMatrix_c A );
ElError ElConjugate_z( ElMatrix_z A );

ElError ElConjugateDist_c( ElDistMatrix_c A );
ElError ElConjugateDist_z( ElDistMatrix_z A );

/* B = A 
   ===== */
ElError ElCopy_i( ElConstMatrix_i A, ElMatrix_i B );
ElError ElCopy_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElCopy_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElCopy_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElCopy_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElCopyDist_i( ElConstDistMatrix_i A, ElDistMatrix_i B );
ElError ElCopyDist_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElCopyDist_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElCopyDist_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElCopyDist_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* DiagonalScale 
   ============= */
ElError ElDiagonalScale_i
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_i d, ElMatrix_i X );
ElError ElDiagonalScale_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X );
ElError ElDiagonalScale_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X );
ElError ElDiagonalScale_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X );
ElError ElDiagonalScale_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X );

ElError ElDiagonalScaleDist_i
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_i d, ElDistMatrix_i X );
ElError ElDiagonalScaleDist_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X );
ElError ElDiagonalScaleDist_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X );
ElError ElDiagonalScaleDist_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X );
ElError ElDiagonalScaleDist_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X );

/* DiagonalScaleTrapezoid
   ====================== */
ElError ElDiagonalScaleTrapezoid_i
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_i d, ElMatrix_i X, ElInt offset );
ElError ElDiagonalScaleTrapezoid_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X, ElInt offset );
ElError ElDiagonalScaleTrapezoid_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X, ElInt offset );
ElError ElDiagonalScaleTrapezoid_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X, ElInt offset );
ElError ElDiagonalScaleTrapezoid_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X, ElInt offset );

ElError ElDiagonalScaleTrapezoidDist_i
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_i d, ElDistMatrix_i X, ElInt offset );
ElError ElDiagonalScaleTrapezoidDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X, ElInt offset );
ElError ElDiagonalScaleTrapezoidDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X, ElInt offset );
ElError ElDiagonalScaleTrapezoidDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X, ElInt offset );
ElError ElDiagonalScaleTrapezoidDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X, ElInt offset );

/* DiagonalSolve 
   ============= */
ElError ElDiagonalSolve_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X );
ElError ElDiagonalSolve_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X );
ElError ElDiagonalSolve_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X );
ElError ElDiagonalSolve_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X );

ElError ElDiagonalSolveDist_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X );
ElError ElDiagonalSolveDist_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X );
ElError ElDiagonalSolveDist_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X );
ElError ElDiagonalSolveDist_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X );

/* DiagonalSolveTrapezoid
   ====================== */
ElError ElDiagonalSolveTrapezoid_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X, ElInt offset );
ElError ElDiagonalSolveTrapezoid_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X, ElInt offset );
ElError ElDiagonalSolveTrapezoid_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X, ElInt offset );
ElError ElDiagonalSolveTrapezoid_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X, ElInt offset );

ElError ElDiagonalSolveTrapezoidDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X, ElInt offset );
ElError ElDiagonalSolveTrapezoidDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X, ElInt offset );
ElError ElDiagonalSolveTrapezoidDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X, ElInt offset );
ElError ElDiagonalSolveTrapezoidDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X, ElInt offset );

/* Dot 
   === */
ElError ElDot_i( ElConstMatrix_i A, ElConstMatrix_i B, ElInt* prod );
ElError ElDot_s( ElConstMatrix_s A, ElConstMatrix_s B, float* prod );
ElError ElDot_d( ElConstMatrix_d A, ElConstMatrix_d B, double* prod );
ElError ElDot_c( ElConstMatrix_c A, ElConstMatrix_c B, complex_float* prod );
ElError ElDot_z( ElConstMatrix_z A, ElConstMatrix_z B, complex_double* prod );

ElError ElDotDist_i
( ElConstDistMatrix_i A, ElConstDistMatrix_i B, ElInt* prod );
ElError ElDotDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, float* prod );
ElError ElDotDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, double* prod );
ElError ElDotDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, complex_float* prod );
ElError ElDotDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, complex_double* prod );

/* Dotu
   ==== */
ElError ElDotu_i( ElConstMatrix_i A, ElConstMatrix_i B, ElInt* prod );
ElError ElDotu_s( ElConstMatrix_s A, ElConstMatrix_s B, float* prod );
ElError ElDotu_d( ElConstMatrix_d A, ElConstMatrix_d B, double* prod );
ElError ElDotu_c( ElConstMatrix_c A, ElConstMatrix_c B, complex_float* prod );
ElError ElDotu_z( ElConstMatrix_z A, ElConstMatrix_z B, complex_double* prod );

ElError ElDotuDist_i
( ElConstDistMatrix_i A, ElConstDistMatrix_i B, ElInt* prod );
ElError ElDotuDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, float* prod );
ElError ElDotuDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, double* prod );
ElError ElDotuDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, complex_float* prod );
ElError ElDotuDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, complex_double* prod );

/* EntrywiseFill
   ============= */
ElError ElEntrywiseFill_i( ElMatrix_i A, ElInt (*fill)() );
ElError ElEntrywiseFill_s( ElMatrix_s A, float (*fill)() );
ElError ElEntrywiseFill_d( ElMatrix_d A, double (*fill)() );
ElError ElEntrywiseFill_c( ElMatrix_c A, complex_float (*fill)() );
ElError ElEntrywiseFill_z( ElMatrix_z A, complex_double (*fill)() );

ElError ElEntrywiseFillDist_i( ElDistMatrix_i A, ElInt (*fill)() );
ElError ElEntrywiseFillDist_s( ElDistMatrix_s A, float (*fill)() );
ElError ElEntrywiseFillDist_d( ElDistMatrix_d A, double (*fill)() );
ElError ElEntrywiseFillDist_c( ElDistMatrix_c A, complex_float (*fill)() );
ElError ElEntrywiseFillDist_z( ElDistMatrix_z A, complex_double (*fill)() );

/* EntrywiseMap
   ============ */
ElError ElEntrywiseMap_i
( ElMatrix_i A, ElInt (*func)(ElInt) );
ElError ElEntrywiseMap_s
( ElMatrix_s A, float (*func)(float) );
ElError ElEntrywiseMap_d
( ElMatrix_d A, double (*func)(double) );
ElError ElEntrywiseMap_c
( ElMatrix_c A, complex_float (*func)(complex_float) );
ElError ElEntrywiseMap_z
( ElMatrix_z A, complex_double (*func)(complex_double) );

ElError ElEntrywiseMapDist_i
( ElDistMatrix_i A, ElInt (*func)(ElInt) );
ElError ElEntrywiseMapDist_s
( ElDistMatrix_s A, float (*func)(float) );
ElError ElEntrywiseMapDist_d
( ElDistMatrix_d A, double (*func)(double) );
ElError ElEntrywiseMapDist_c
( ElDistMatrix_c A, complex_float (*func)(complex_float) );
ElError ElEntrywiseMapDist_z
( ElDistMatrix_z A, complex_double (*func)(complex_double) );

/* Fill
   ==== */
ElError ElFill_i( ElMatrix_i A, ElInt alpha );
ElError ElFill_s( ElMatrix_s A, float alpha );
ElError ElFill_d( ElMatrix_d A, double alpha );
ElError ElFill_c( ElMatrix_c A, complex_float alpha );
ElError ElFill_z( ElMatrix_z A, complex_double alpha );

ElError ElFillDist_i( ElDistMatrix_i A, ElInt alpha );
ElError ElFillDist_s( ElDistMatrix_s A, float alpha );
ElError ElFillDist_d( ElDistMatrix_d A, double alpha );
ElError ElFillDist_c( ElDistMatrix_c A, complex_float alpha );
ElError ElFillDist_z( ElDistMatrix_z A, complex_double alpha );

/* Hadamard
   ======== */
ElError ElHadamard_i( ElConstMatrix_i A, ElConstMatrix_i B, ElMatrix_i C );
ElError ElHadamard_s( ElConstMatrix_s A, ElConstMatrix_s B, ElMatrix_s C );
ElError ElHadamard_d( ElConstMatrix_d A, ElConstMatrix_d B, ElMatrix_d C );
ElError ElHadamard_c( ElConstMatrix_c A, ElConstMatrix_c B, ElMatrix_c C );
ElError ElHadamard_z( ElConstMatrix_z A, ElConstMatrix_z B, ElMatrix_z C );

ElError ElHadamardDist_i
( ElConstDistMatrix_i A, ElConstDistMatrix_i B, ElDistMatrix_i C );
ElError ElHadamardDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, ElDistMatrix_s C );
ElError ElHadamardDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, ElDistMatrix_d C );
ElError ElHadamardDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, ElDistMatrix_c C );
ElError ElHadamardDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, ElDistMatrix_z C );

/* HilbertSchmidt
   ============== */
/* NOTE: This is the same as Dot */
ElError ElHilbertSchmidt_i
( ElConstMatrix_i A, ElConstMatrix_i B, ElInt* prod );
ElError ElHilbertSchmidt_s
( ElConstMatrix_s A, ElConstMatrix_s B, float* prod );
ElError ElHilbertSchmidt_d
( ElConstMatrix_d A, ElConstMatrix_d B, double* prod );
ElError ElHilbertSchmidt_c
( ElConstMatrix_c A, ElConstMatrix_c B, complex_float* prod );
ElError ElHilbertSchmidt_z
( ElConstMatrix_z A, ElConstMatrix_z B, complex_double* prod );

ElError ElHilbertSchmidtDist_i
( ElConstDistMatrix_i A, ElConstDistMatrix_i B, ElInt* prod );
ElError ElHilbertSchmidtDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, float* prod );
ElError ElHilbertSchmidtDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, double* prod );
ElError ElHilbertSchmidtDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, complex_float* prod );
ElError ElHilbertSchmidtDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, complex_double* prod );

/* IndexDependentFill
   ================== */
ElError ElIndexDependentFill_i
( ElMatrix_i A, ElInt (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFill_s
( ElMatrix_s A, float (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFill_d
( ElMatrix_d A, double (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFill_c
( ElMatrix_c A, complex_float (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFill_z
( ElMatrix_z A, complex_double (*fill)(ElInt,ElInt) );

ElError ElIndexDependentFillDist_i
( ElDistMatrix_i A, ElInt (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillDist_s
( ElDistMatrix_s A, float (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillDist_d
( ElDistMatrix_d A, double (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillDist_c
( ElDistMatrix_c A, complex_float (*fill)(ElInt,ElInt) );
ElError ElIndexDependentFillDist_z
( ElDistMatrix_z A, complex_double (*fill)(ElInt,ElInt) );

/* IndexDependentMap
   ================= */
ElError ElIndexDependentMap_i
( ElMatrix_i A, ElInt (*func)(ElInt,ElInt,ElInt) );
ElError ElIndexDependentMap_s
( ElMatrix_s A, float (*func)(ElInt,ElInt,float) );
ElError ElIndexDependentMap_d
( ElMatrix_d A, double (*func)(ElInt,ElInt,double) );
ElError ElIndexDependentMap_c
( ElMatrix_c A, complex_float (*func)(ElInt,ElInt,complex_float) );
ElError ElIndexDependentMap_z
( ElMatrix_z A, complex_double (*func)(ElInt,ElInt,complex_double) );

ElError ElIndexDependentMapDist_i
( ElDistMatrix_i A, ElInt (*func)(ElInt,ElInt,ElInt) );
ElError ElIndexDependentMapDist_s
( ElDistMatrix_s A, float (*func)(ElInt,ElInt,float) );
ElError ElIndexDependentMapDist_d
( ElDistMatrix_d A, double (*func)(ElInt,ElInt,double) );
ElError ElIndexDependentMapDist_c
( ElDistMatrix_c A, complex_float (*func)(ElInt,ElInt,complex_float) );
ElError ElIndexDependentMapDist_z
( ElDistMatrix_z A, complex_double (*func)(ElInt,ElInt,complex_double) );

/* MakeHermitian 
   ============= */
ElError ElMakeHermitian_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElMakeHermitian_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElMakeHermitianDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElMakeHermitianDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* MakeReal
   ======== */
ElError ElMakeReal_c( ElMatrix_c A );
ElError ElMakeReal_z( ElMatrix_z A );

ElError ElMakeRealDist_c( ElDistMatrix_c A );
ElError ElMakeRealDist_z( ElDistMatrix_z A );

/* MakeSymmetric 
   ============= */
ElError ElMakeSymmetric_i( ElUpperOrLower uplo, ElMatrix_i A );
ElError ElMakeSymmetric_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElMakeSymmetric_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElMakeSymmetric_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElMakeSymmetric_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElMakeSymmetricDist_i( ElUpperOrLower uplo, ElDistMatrix_i A );
ElError ElMakeSymmetricDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElMakeSymmetricDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElMakeSymmetricDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElMakeSymmetricDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* MakeTrapezoidal
   =============== */
ElError ElMakeTrapezoidal_i( ElUpperOrLower uplo, ElMatrix_i A, ElInt offset );
ElError ElMakeTrapezoidal_s( ElUpperOrLower uplo, ElMatrix_s A, ElInt offset );
ElError ElMakeTrapezoidal_d( ElUpperOrLower uplo, ElMatrix_d A, ElInt offset );
ElError ElMakeTrapezoidal_c( ElUpperOrLower uplo, ElMatrix_c A, ElInt offset );
ElError ElMakeTrapezoidal_z( ElUpperOrLower uplo, ElMatrix_z A, ElInt offset );

ElError ElMakeTrapezoidalDist_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElInt offset );
ElError ElMakeTrapezoidalDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElInt offset );
ElError ElMakeTrapezoidalDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElInt offset );
ElError ElMakeTrapezoidalDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElInt offset );
ElError ElMakeTrapezoidalDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElInt offset );

/* MakeTriangular
   ============== */
ElError ElMakeTriangular_i( ElUpperOrLower uplo, ElMatrix_i A );
ElError ElMakeTriangular_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElMakeTriangular_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElMakeTriangular_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElMakeTriangular_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElMakeTriangularDist_i( ElUpperOrLower uplo, ElDistMatrix_i A );
ElError ElMakeTriangularDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElMakeTriangularDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElMakeTriangularDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElMakeTriangularDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Max
   === */
ElError ElMax_i( ElConstMatrix_i A, ElValueIntPair_i* entry );
ElError ElMax_s( ElConstMatrix_s A, ElValueIntPair_s* entry );
ElError ElMax_d( ElConstMatrix_d A, ElValueIntPair_d* entry );

ElError ElMaxDist_i( ElConstDistMatrix_i A, ElValueIntPair_i* entry );
ElError ElMaxDist_s( ElConstDistMatrix_s A, ElValueIntPair_s* entry );
ElError ElMaxDist_d( ElConstDistMatrix_d A, ElValueIntPair_d* entry );

ElError ElSymmetricMax_i
( ElUpperOrLower uplo, ElConstMatrix_i A, ElValueIntPair_i* entry );
ElError ElSymmetricMax_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElValueIntPair_s* entry );
ElError ElSymmetricMax_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElValueIntPair_d* entry );

ElError ElSymmetricMaxDist_i
( ElUpperOrLower uplo, ElConstDistMatrix_i A, ElValueIntPair_i* entry );
ElError ElSymmetricMaxDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElValueIntPair_s* entry );
ElError ElSymmetricMaxDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElValueIntPair_d* entry );

ElError ElVectorMax_i( ElConstMatrix_i x, ElValueInt_i* entry );
ElError ElVectorMax_s( ElConstMatrix_s x, ElValueInt_s* entry );
ElError ElVectorMax_d( ElConstMatrix_d x, ElValueInt_d* entry );

ElError ElVectorMaxDist_i( ElConstDistMatrix_i x, ElValueInt_i* entry );
ElError ElVectorMaxDist_s( ElConstDistMatrix_s x, ElValueInt_s* entry );
ElError ElVectorMaxDist_d( ElConstDistMatrix_d x, ElValueInt_d* entry );

/* MaxAbs
   ====== */
ElError ElMaxAbs_i( ElConstMatrix_i A, ElValueIntPair_i* entry );
ElError ElMaxAbs_s( ElConstMatrix_s A, ElValueIntPair_s* entry );
ElError ElMaxAbs_d( ElConstMatrix_d A, ElValueIntPair_d* entry );
ElError ElMaxAbs_c( ElConstMatrix_c A, ElValueIntPair_s* entry );
ElError ElMaxAbs_z( ElConstMatrix_z A, ElValueIntPair_d* entry );

ElError ElMaxAbsDist_i( ElConstDistMatrix_i A, ElValueIntPair_i* entry );
ElError ElMaxAbsDist_s( ElConstDistMatrix_s A, ElValueIntPair_s* entry );
ElError ElMaxAbsDist_d( ElConstDistMatrix_d A, ElValueIntPair_d* entry );
ElError ElMaxAbsDist_c( ElConstDistMatrix_c A, ElValueIntPair_s* entry );
ElError ElMaxAbsDist_z( ElConstDistMatrix_z A, ElValueIntPair_d* entry );

ElError ElSymmetricMaxAbs_i
( ElUpperOrLower uplo, ElConstMatrix_i A, ElValueIntPair_i* entry );
ElError ElSymmetricMaxAbs_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElValueIntPair_s* entry );
ElError ElSymmetricMaxAbs_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElValueIntPair_d* entry );
ElError ElSymmetricMaxAbs_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElValueIntPair_s* entry );
ElError ElSymmetricMaxAbs_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElValueIntPair_d* entry );

ElError ElSymmetricMaxAbsDist_i
( ElUpperOrLower uplo, ElConstDistMatrix_i A, ElValueIntPair_i* entry );
ElError ElSymmetricMaxAbsDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElValueIntPair_s* entry );
ElError ElSymmetricMaxAbsDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElValueIntPair_d* entry );
ElError ElSymmetricMaxAbsDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElValueIntPair_s* entry );
ElError ElSymmetricMaxAbsDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElValueIntPair_d* entry );

ElError ElVectorMaxAbs_i( ElConstMatrix_i x, ElValueInt_i* entry );
ElError ElVectorMaxAbs_s( ElConstMatrix_s x, ElValueInt_s* entry );
ElError ElVectorMaxAbs_d( ElConstMatrix_d x, ElValueInt_d* entry );
ElError ElVectorMaxAbs_c( ElConstMatrix_c x, ElValueInt_s* entry );
ElError ElVectorMaxAbs_z( ElConstMatrix_z x, ElValueInt_d* entry );

ElError ElVectorMaxAbsDist_i( ElConstDistMatrix_i x, ElValueInt_i* entry );
ElError ElVectorMaxAbsDist_s( ElConstDistMatrix_s x, ElValueInt_s* entry );
ElError ElVectorMaxAbsDist_d( ElConstDistMatrix_d x, ElValueInt_d* entry );
ElError ElVectorMaxAbsDist_c( ElConstDistMatrix_c x, ElValueInt_s* entry );
ElError ElVectorMaxAbsDist_z( ElConstDistMatrix_z x, ElValueInt_d* entry );

/* Min
   === */
ElError ElMin_i( ElConstMatrix_i A, ElValueIntPair_i* entry );
ElError ElMin_s( ElConstMatrix_s A, ElValueIntPair_s* entry );
ElError ElMin_d( ElConstMatrix_d A, ElValueIntPair_d* entry );

ElError ElMinDist_i( ElConstDistMatrix_i A, ElValueIntPair_i* entry );
ElError ElMinDist_s( ElConstDistMatrix_s A, ElValueIntPair_s* entry );
ElError ElMinDist_d( ElConstDistMatrix_d A, ElValueIntPair_d* entry );

ElError ElSymmetricMin_i
( ElUpperOrLower uplo, ElConstMatrix_i A, ElValueIntPair_i* entry );
ElError ElSymmetricMin_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElValueIntPair_s* entry );
ElError ElSymmetricMin_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElValueIntPair_d* entry );

ElError ElSymmetricMinDist_i
( ElUpperOrLower uplo, ElConstDistMatrix_i A, ElValueIntPair_i* entry );
ElError ElSymmetricMinDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElValueIntPair_s* entry );
ElError ElSymmetricMinDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElValueIntPair_d* entry );

ElError ElVectorMin_i( ElConstMatrix_i x, ElValueInt_i* entry );
ElError ElVectorMin_s( ElConstMatrix_s x, ElValueInt_s* entry );
ElError ElVectorMin_d( ElConstMatrix_d x, ElValueInt_d* entry );

ElError ElVectorMinDist_i( ElConstDistMatrix_i x, ElValueInt_i* entry );
ElError ElVectorMinDist_s( ElConstDistMatrix_s x, ElValueInt_s* entry );
ElError ElVectorMinDist_d( ElConstDistMatrix_d x, ElValueInt_d* entry );

/* MinAbs
   ====== */
ElError ElMinAbs_i( ElConstMatrix_i A, ElValueIntPair_i* entry );
ElError ElMinAbs_s( ElConstMatrix_s A, ElValueIntPair_s* entry );
ElError ElMinAbs_d( ElConstMatrix_d A, ElValueIntPair_d* entry );
ElError ElMinAbs_c( ElConstMatrix_c A, ElValueIntPair_s* entry );
ElError ElMinAbs_z( ElConstMatrix_z A, ElValueIntPair_d* entry );

ElError ElMinAbsDist_i( ElConstDistMatrix_i A, ElValueIntPair_i* entry );
ElError ElMinAbsDist_s( ElConstDistMatrix_s A, ElValueIntPair_s* entry );
ElError ElMinAbsDist_d( ElConstDistMatrix_d A, ElValueIntPair_d* entry );
ElError ElMinAbsDist_c( ElConstDistMatrix_c A, ElValueIntPair_s* entry );
ElError ElMinAbsDist_z( ElConstDistMatrix_z A, ElValueIntPair_d* entry );

ElError ElSymmetricMinAbs_i
( ElUpperOrLower uplo, ElConstMatrix_i A, ElValueIntPair_i* entry );
ElError ElSymmetricMinAbs_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElValueIntPair_s* entry );
ElError ElSymmetricMinAbs_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElValueIntPair_d* entry );
ElError ElSymmetricMinAbs_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElValueIntPair_s* entry );
ElError ElSymmetricMinAbs_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElValueIntPair_d* entry );

ElError ElSymmetricMinAbsDist_i
( ElUpperOrLower uplo, ElConstDistMatrix_i A, ElValueIntPair_i* entry );
ElError ElSymmetricMinAbsDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElValueIntPair_s* entry );
ElError ElSymmetricMinAbsDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElValueIntPair_d* entry );
ElError ElSymmetricMinAbsDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElValueIntPair_s* entry );
ElError ElSymmetricMinAbsDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElValueIntPair_d* entry );

ElError ElVectorMinAbs_i( ElConstMatrix_i x, ElValueInt_i* entry );
ElError ElVectorMinAbs_s( ElConstMatrix_s x, ElValueInt_s* entry );
ElError ElVectorMinAbs_d( ElConstMatrix_d x, ElValueInt_d* entry );
ElError ElVectorMinAbs_c( ElConstMatrix_c x, ElValueInt_s* entry );
ElError ElVectorMinAbs_z( ElConstMatrix_z x, ElValueInt_d* entry );

ElError ElVectorMinAbsDist_i( ElConstDistMatrix_i x, ElValueInt_i* entry );
ElError ElVectorMinAbsDist_s( ElConstDistMatrix_s x, ElValueInt_s* entry );
ElError ElVectorMinAbsDist_d( ElConstDistMatrix_d x, ElValueInt_d* entry );
ElError ElVectorMinAbsDist_c( ElConstDistMatrix_c x, ElValueInt_s* entry );
ElError ElVectorMinAbsDist_z( ElConstDistMatrix_z x, ElValueInt_d* entry );

/* Nrm2
   ==== */
ElError ElNrm2_s( ElConstMatrix_s A, float* gamma );
ElError ElNrm2_d( ElConstMatrix_d A, double* gamma );
ElError ElNrm2_c( ElConstMatrix_c A, float* gamma );
ElError ElNrm2_z( ElConstMatrix_z A, double* gamma );

ElError ElNrm2Dist_s( ElConstDistMatrix_s A, float* gamma );
ElError ElNrm2Dist_d( ElConstDistMatrix_d A, double* gamma );
ElError ElNrm2Dist_c( ElConstDistMatrix_c A, float* gamma );
ElError ElNrm2Dist_z( ElConstDistMatrix_z A, double* gamma );

/* Scale
   ===== */
ElError ElScale_i( ElInt alpha, ElMatrix_i A );
ElError ElScale_s( float alpha, ElMatrix_s A );
ElError ElScale_d( double alpha, ElMatrix_d A );
ElError ElScale_c( complex_float alpha, ElMatrix_c A );
ElError ElScale_z( complex_double alpha, ElMatrix_z A );

ElError ElScaleDist_i( ElInt alpha, ElDistMatrix_i A );
ElError ElScaleDist_s( float alpha, ElDistMatrix_s A );
ElError ElScaleDist_d( double alpha, ElDistMatrix_d A );
ElError ElScaleDist_c( complex_float alpha, ElDistMatrix_c A );
ElError ElScaleDist_z( complex_double alpha, ElDistMatrix_z A );

/* ScaleTrapezoid
   ============== */
ElError ElScaleTrapezoid_i
( ElInt alpha, ElUpperOrLower uplo, ElMatrix_i A, ElInt offset );
ElError ElScaleTrapezoid_s
( float alpha, ElUpperOrLower uplo, ElMatrix_s A, ElInt offset );
ElError ElScaleTrapezoid_d
( double alpha, ElUpperOrLower uplo, ElMatrix_d A, ElInt offset );
ElError ElScaleTrapezoid_c
( complex_float alpha, ElUpperOrLower uplo, ElMatrix_c A, ElInt offset );
ElError ElScaleTrapezoid_z
( complex_double alpha, ElUpperOrLower uplo, ElMatrix_z A, ElInt offset );

ElError ElScaleTrapezoidDist_i
( ElInt alpha, ElUpperOrLower uplo, ElDistMatrix_i A, ElInt offset );
ElError ElScaleTrapezoidDist_s
( float alpha, ElUpperOrLower uplo, ElDistMatrix_s A, ElInt offset );
ElError ElScaleTrapezoidDist_d
( double alpha, ElUpperOrLower uplo, ElDistMatrix_d A, ElInt offset );
ElError ElScaleTrapezoidDist_c
( complex_float alpha, ElUpperOrLower uplo, ElDistMatrix_c A, ElInt offset );
ElError ElScaleTrapezoidDist_z
( complex_double alpha, ElUpperOrLower uplo, ElDistMatrix_z A, ElInt offset );

/* SetDiagonal
   =========== */
ElError ElSetDiagonal_i( ElMatrix_i A, ElInt alpha, ElInt offset );
ElError ElSetDiagonal_s( ElMatrix_s A, float alpha, ElInt offset );
ElError ElSetDiagonal_d( ElMatrix_d A, double alpha, ElInt offset );
ElError ElSetDiagonal_c( ElMatrix_c A, complex_float alpha, ElInt offset );
ElError ElSetDiagonal_z( ElMatrix_z A, complex_double alpha, ElInt offset );

ElError ElSetDiagonalDist_i
( ElDistMatrix_i A, ElInt alpha, ElInt offset );
ElError ElSetDiagonalDist_s
( ElDistMatrix_s A, float alpha, ElInt offset );
ElError ElSetDiagonalDist_d
( ElDistMatrix_d A, double alpha, ElInt offset );
ElError ElSetDiagonalDist_c
( ElDistMatrix_c A, complex_float alpha, ElInt offset );
ElError ElSetDiagonalDist_z
( ElDistMatrix_z A, complex_double alpha, ElInt offset );

/* Swap
   ==== */
ElError ElSwap_i( ElOrientation orientation, ElMatrix_i X, ElMatrix_i Y );
ElError ElSwap_s( ElOrientation orientation, ElMatrix_s X, ElMatrix_s Y );
ElError ElSwap_d( ElOrientation orientation, ElMatrix_d X, ElMatrix_d Y );
ElError ElSwap_c( ElOrientation orientation, ElMatrix_c X, ElMatrix_c Y );
ElError ElSwap_z( ElOrientation orientation, ElMatrix_z X, ElMatrix_z Y );

ElError ElSwapDist_i
( ElOrientation orientation, ElDistMatrix_i X, ElDistMatrix_i Y );
ElError ElSwapDist_s
( ElOrientation orientation, ElDistMatrix_s X, ElDistMatrix_s Y );
ElError ElSwapDist_d
( ElOrientation orientation, ElDistMatrix_d X, ElDistMatrix_d Y );
ElError ElSwapDist_c
( ElOrientation orientation, ElDistMatrix_c X, ElDistMatrix_c Y );
ElError ElSwapDist_z
( ElOrientation orientation, ElDistMatrix_z X, ElDistMatrix_z Y );

ElError ElRowSwap_i( ElMatrix_i A, ElInt to, ElInt from );
ElError ElRowSwap_s( ElMatrix_s A, ElInt to, ElInt from );
ElError ElRowSwap_d( ElMatrix_d A, ElInt to, ElInt from );
ElError ElRowSwap_c( ElMatrix_c A, ElInt to, ElInt from );
ElError ElRowSwap_z( ElMatrix_z A, ElInt to, ElInt from );

ElError ElRowSwapDist_i( ElDistMatrix_i A, ElInt to, ElInt from );
ElError ElRowSwapDist_s( ElDistMatrix_s A, ElInt to, ElInt from );
ElError ElRowSwapDist_d( ElDistMatrix_d A, ElInt to, ElInt from );
ElError ElRowSwapDist_c( ElDistMatrix_c A, ElInt to, ElInt from );
ElError ElRowSwapDist_z( ElDistMatrix_z A, ElInt to, ElInt from );

ElError ElColSwap_i( ElMatrix_i A, ElInt to, ElInt from );
ElError ElColSwap_s( ElMatrix_s A, ElInt to, ElInt from );
ElError ElColSwap_d( ElMatrix_d A, ElInt to, ElInt from );
ElError ElColSwap_c( ElMatrix_c A, ElInt to, ElInt from );
ElError ElColSwap_z( ElMatrix_z A, ElInt to, ElInt from );

ElError ElColSwapDist_i( ElDistMatrix_i A, ElInt to, ElInt from );
ElError ElColSwapDist_s( ElDistMatrix_s A, ElInt to, ElInt from );
ElError ElColSwapDist_d( ElDistMatrix_d A, ElInt to, ElInt from );
ElError ElColSwapDist_c( ElDistMatrix_c A, ElInt to, ElInt from );
ElError ElColSwapDist_z( ElDistMatrix_z A, ElInt to, ElInt from );

ElError ElSymmetricSwap_i
( ElUpperOrLower uplo, ElMatrix_i A, ElInt to, ElInt from );
ElError ElSymmetricSwap_s
( ElUpperOrLower uplo, ElMatrix_s A, ElInt to, ElInt from );
ElError ElSymmetricSwap_d
( ElUpperOrLower uplo, ElMatrix_d A, ElInt to, ElInt from );
ElError ElSymmetricSwap_c
( ElUpperOrLower uplo, ElMatrix_c A, ElInt to, ElInt from );
ElError ElSymmetricSwap_z
( ElUpperOrLower uplo, ElMatrix_z A, ElInt to, ElInt from );

ElError ElSymmetricSwapDist_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElInt to, ElInt from );
ElError ElSymmetricSwapDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElInt to, ElInt from );
ElError ElSymmetricSwapDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElInt to, ElInt from );
ElError ElSymmetricSwapDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElInt to, ElInt from );
ElError ElSymmetricSwapDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElInt to, ElInt from );

ElError ElHermitianSwap_c
( ElUpperOrLower uplo, ElMatrix_c A, ElInt to, ElInt from );
ElError ElHermitianSwap_z
( ElUpperOrLower uplo, ElMatrix_z A, ElInt to, ElInt from );

ElError ElHermitianSwapDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElInt to, ElInt from );
ElError ElHermitianSwapDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElInt to, ElInt from );

/* B = A^T 
   ======= */
ElError ElTranspose_i( ElConstMatrix_i A, ElMatrix_i B );
ElError ElTranspose_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElTranspose_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElTranspose_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElTranspose_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElTransposeDist_i( ElConstDistMatrix_i A, ElDistMatrix_i B );
ElError ElTransposeDist_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElTransposeDist_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElTransposeDist_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElTransposeDist_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* UpdateDiagonal
   ============== */
ElError ElUpdateDiagonal_i( ElMatrix_i A, ElInt alpha, ElInt offset );
ElError ElUpdateDiagonal_s( ElMatrix_s A, float alpha, ElInt offset );
ElError ElUpdateDiagonal_d( ElMatrix_d A, double alpha, ElInt offset );
ElError ElUpdateDiagonal_c( ElMatrix_c A, complex_float alpha, ElInt offset );
ElError ElUpdateDiagonal_z( ElMatrix_z A, complex_double alpha, ElInt offset );

ElError ElUpdateDiagonalDist_i
( ElDistMatrix_i A, ElInt alpha, ElInt offset );
ElError ElUpdateDiagonalDist_s
( ElDistMatrix_s A, float alpha, ElInt offset );
ElError ElUpdateDiagonalDist_d
( ElDistMatrix_d A, double alpha, ElInt offset );
ElError ElUpdateDiagonalDist_c
( ElDistMatrix_c A, complex_float alpha, ElInt offset );
ElError ElUpdateDiagonalDist_z
( ElDistMatrix_z A, complex_double alpha, ElInt offset );

/* Zero
   ==== */
ElError ElZero_i( ElMatrix_i A );
ElError ElZero_s( ElMatrix_s A );
ElError ElZero_d( ElMatrix_d A );
ElError ElZero_c( ElMatrix_c A );
ElError ElZero_z( ElMatrix_z A );

ElError ElZeroDist_i( ElDistMatrix_i A );
ElError ElZeroDist_s( ElDistMatrix_s A );
ElError ElZeroDist_d( ElDistMatrix_d A );
ElError ElZeroDist_c( ElDistMatrix_c A );
ElError ElZeroDist_z( ElDistMatrix_z A );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_BLAS_LEVEL1_C_H */
