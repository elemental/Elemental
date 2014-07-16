/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El-C.h"
using namespace El;

extern "C" {

#define C_PROTO_BASE(SIG,T) \
  /* B = A^H */ \
  ElError ElAdjointMatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Adjoint( *Reinterpret(A), *Reinterpret(B) ) ) } \
  ElError ElAdjointDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Adjoint( *Reinterpret(A), *Reinterpret(B) ) ) } \
  /* Y := alpha X + Y */ \
  ElError ElAxpyMatrix_ ## SIG \
  ( CREFLECT(T) alpha, ElConstMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY( Axpy( Reinterpret(alpha), *Reinterpret(X), *Reinterpret(Y) ) ) } \
  ElError ElAxpyDistMatrix_ ## SIG \
  ( CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY( Axpy( Reinterpret(alpha), *Reinterpret(X), *Reinterpret(Y) ) ) } \
  /* tri(Y) := tri(alpha X + Y) */ \
  ElError ElAxpyTriangleMatrix_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY \
    ( AxpyTriangle \
      ( Reinterpret(uplo), Reinterpret(alpha), \
        *Reinterpret(X), *Reinterpret(Y) ) ) } \
  ElError ElAxpyTriangleDistMatrix_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY \
    ( AxpyTriangle \
      ( Reinterpret(uplo), Reinterpret(alpha), \
        *Reinterpret(X), *Reinterpret(Y) ) ) } \
  /* B = A */ \
  ElError ElCopyMatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *Reinterpret(A), *Reinterpret(B) ) ) } \
  ElError ElCopyDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *Reinterpret(A), *Reinterpret(B) ) ) } \
  /* DiagonalScale */ \
  ElError ElDiagonalScaleMatrix_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG X ) \
  { EL_TRY( \
      DiagonalScale \
      ( Reinterpret(side), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X) ) ) } \
  ElError ElDiagonalScaleDistMatrix_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( \
      DiagonalScale \
      ( Reinterpret(side), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X) ) ) } \
  /* DiagonalScaleTrapezoid */ \
  ElError ElDiagonalScaleTrapezoidMatrix_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG X, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( Reinterpret(side), Reinterpret(uplo), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X), offset ) ) } \
  ElError ElDiagonalScaleTrapezoidDistMatrix_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG X, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( Reinterpret(side), Reinterpret(uplo), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X), offset ) ) } \
  /* Dot product (<A,B>=vec(A)^H vec(B)) */ \
  ElError ElDotMatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = Reinterpret(Dot(*Reinterpret(A),*Reinterpret(B))) ) } \
  ElError ElDotDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = Reinterpret(Dot(*Reinterpret(A),*Reinterpret(B))) ) } \
  /* Unconjugated dot product, vec(A)^T vec(B) */ \
  ElError ElDotuMatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = Reinterpret(Dotu(*Reinterpret(A),*Reinterpret(B))) ) } \
  ElError ElDotuDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = Reinterpret(Dotu(*Reinterpret(A),*Reinterpret(B))) ) } \
  /* EntrywiseFill */ \
  ElError ElEntrywiseFillMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*fill)() ) \
  { try { \
      auto newFill = [&]() { return Reinterpret(fill()); }; \
      EntrywiseFill( *Reinterpret(A), std::function<T(void)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseFillDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*fill)() ) \
  { try { \
      auto newFill = [&]() { return Reinterpret(fill()); }; \
      EntrywiseFill( *Reinterpret(A), std::function<T(void)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* EntrywiseMap */ \
  ElError ElEntrywiseMapMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return Reinterpret(func(Reinterpret(alpha))); }; \
      EntrywiseMap( *Reinterpret(A), std::function<T(T)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseMapDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return Reinterpret(func(Reinterpret(alpha))); }; \
      EntrywiseMap( *Reinterpret(A), std::function<T(T)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Fill */ \
  ElError ElFillMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *Reinterpret(A), Reinterpret(alpha) ) ) } \
  ElError ElFillDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *Reinterpret(A), Reinterpret(alpha) ) ) } \
  /* Hadamard */ \
  ElError ElHadamardMatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG C ) \
  { EL_TRY( Hadamard(*Reinterpret(A),*Reinterpret(B),*Reinterpret(C)) ) } \
  ElError ElHadamardDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG C ) \
  { EL_TRY( Hadamard(*Reinterpret(A),*Reinterpret(B),*Reinterpret(C)) ) } \
  /* Hilbert-Schmidt inner product (same as Dot) */ \
  ElError ElHilbertSchmidtMatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = \
      Reinterpret(HilbertSchmidt(*Reinterpret(A),*Reinterpret(B))) ) } \
  ElError ElHilbertSchmidtDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = \
      Reinterpret(HilbertSchmidt(*Reinterpret(A),*Reinterpret(B))) ) } \
  /* IndexDependentFill */ \
  ElError ElIndexDependentFillMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*fill)(ElInt,ElInt) ) \
  { try { \
      auto newFill = [&]( Int i, Int j ) { return Reinterpret(fill(i,j)); }; \
      IndexDependentFill \
      ( *Reinterpret(A), std::function<T(Int,Int)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElIndexDependentFillDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*fill)(ElInt,ElInt) ) \
  { try { \
      auto newFill = [&]( Int i, Int j ) { return Reinterpret(fill(i,j)); }; \
      IndexDependentFill \
      ( *Reinterpret(A), std::function<T(Int,Int)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* IndexDependentMap */ \
  ElError ElIndexDependentMapMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*func)(ElInt,ElInt,CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( Int i, Int j, T alpha ) \
        { return Reinterpret(func(i,j,Reinterpret(alpha))); }; \
      IndexDependentMap \
      ( *Reinterpret(A), std::function<T(Int,Int,T)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElIndexDependentMapDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*func)(ElInt,ElInt,CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( Int i, Int j, T alpha ) \
        { return Reinterpret(func(i,j,Reinterpret(alpha))); }; \
      IndexDependentMap \
      ( *Reinterpret(A), std::function<T(Int,Int,T)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* MakeHermitian */ \
  ElError ElMakeHermitianMatrix_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( MakeHermitian( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElMakeHermitianDistMatrix_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( MakeHermitian( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  /* MakeSymmetric */ \
  ElError ElMakeSymmetricMatrix_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( MakeSymmetric( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElMakeSymmetricDistMatrix_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( MakeSymmetric( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  /* MakeTrapezoidal */ \
  ElError ElMakeTrapezoidalMatrix_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( MakeTrapezoidal( Reinterpret(uplo), *Reinterpret(A), offset ) ) } \
  ElError ElMakeTrapezoidalDistMatrix_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( MakeTrapezoidal( Reinterpret(uplo), *Reinterpret(A), offset ) ) } \
  /* MakeTriangular */ \
  ElError ElMakeTriangularMatrix_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( MakeTriangular( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElMakeTriangularDistMatrix_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( MakeTriangular( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  /* TODO: Max */ \
  /* TODO: MaxAbs */ \
  /* TODO: Min */ \
  /* TODO: MinAbs */ \
  /* TODO: QuasiDiagonalScale */ \
  /* Scale */ \
  ElError ElScaleMatrix_ ## SIG \
  ( CREFLECT(T) alpha, ElMatrix_ ## SIG A ) \
  { EL_TRY( Scale( Reinterpret(alpha), *Reinterpret(A) ) ) } \
  ElError ElScaleDistMatrix_ ## SIG \
  ( CREFLECT(T) alpha, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Scale( Reinterpret(alpha), *Reinterpret(A) ) ) } \
  /* ScaleTrapezoid */ \
  ElError ElScaleTrapezoidMatrix_ ## SIG \
  ( CREFLECT(T) alpha, ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( \
      ScaleTrapezoid \
      ( Reinterpret(alpha), Reinterpret(uplo), *Reinterpret(A), offset ) ) } \
  ElError ElScaleTrapezoidDistMatrix_ ## SIG \
  ( CREFLECT(T) alpha, ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, \
    ElInt offset ) \
  { EL_TRY( \
      ScaleTrapezoid \
      ( Reinterpret(alpha), Reinterpret(uplo), *Reinterpret(A), offset ) ) } \
  /* SetDiagonal */ \
  ElError ElSetDiagonalMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( SetDiagonal( *Reinterpret(A), Reinterpret(alpha), offset ) ) } \
  ElError ElSetDiagonalDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( SetDiagonal( *Reinterpret(A), Reinterpret(alpha), offset ) ) } \
  /* TODO: Swap */ \
  /* TODO: Symmetric2x2Scale */ \
  /* B = A^T */ \
  ElError ElTransposeMatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*Reinterpret(A),*Reinterpret(B),false) ) } \
  ElError ElTransposeDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*Reinterpret(A),*Reinterpret(B),false) ) } \
  /* UpdateDiagonal */ \
  ElError ElUpdateDiagonalMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( UpdateDiagonal( *Reinterpret(A), Reinterpret(alpha), offset ) ) } \
  ElError ElUpdateDiagonalDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( UpdateDiagonal( *Reinterpret(A), Reinterpret(alpha), offset ) ) } \
  /* Zero */ \
  ElError ElZeroMatrix_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Zero( *Reinterpret(A) ) ) } \
  ElError ElZeroDistMatrix_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Zero( *Reinterpret(A) ) ) } \

#define C_PROTO_INT(SIG,T) C_PROTO_BASE(SIG,T)

#define C_PROTO(SIG,T) \
  C_PROTO_BASE(SIG,T) \
  /* DiagonalSolve */ \
  ElError ElDiagonalSolveMatrix_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG X ) \
  { EL_TRY( \
      DiagonalSolve \
      ( Reinterpret(side), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X) ) ) } \
  ElError ElDiagonalSolveDistMatrix_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( \
      DiagonalSolve \
      ( Reinterpret(side), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X) ) ) } \
  /* Nrm2 (same as FrobeniusNorm) */ \
  ElError ElNrm2Matrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<T> *gamma ) \
  { EL_TRY( *gamma = Nrm2(*Reinterpret(A)) ) } \
  ElError ElNrm2DistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<T> *gamma ) \
  { EL_TRY( *gamma = Nrm2(*Reinterpret(A)) ) } \
  /* TODO: QuasiDiagonalSolve */ \
  /* TODO: Symmetric2x2Inv */ \
  /* TODO: Symmetric2x2Solve */

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO(SIG,T) \
  /* Conjugate */ \
  ElError ElConjugateMatrix_ ## SIG( ElMatrix_ ## SIG A ) \
  { EL_TRY( Conjugate( *Reinterpret(A) ) ) } \
  ElError ElConjugateDistMatrix_ ## SIG( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Conjugate( *Reinterpret(A) ) ) } \
  /* MakeReal */ \
  ElError ElMakeRealMatrix_ ## SIG \
  ( ElMatrix_ ## SIG A ) { EL_TRY( MakeReal( *Reinterpret(A) ) ) } \
  ElError ElMakeRealDistMatrix_ ## SIG \
  ( ElDistMatrix_ ## SIG A ) { EL_TRY( MakeReal( *Reinterpret(A) ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
