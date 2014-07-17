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
  ElError ElAdjoint_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Adjoint( *Reinterpret(A), *Reinterpret(B) ) ) } \
  ElError ElAdjointDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Adjoint( *Reinterpret(A), *Reinterpret(B) ) ) } \
  /* Y := alpha X + Y */ \
  ElError ElAxpy_ ## SIG \
  ( CREFLECT(T) alpha, ElConstMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY( Axpy( Reinterpret(alpha), *Reinterpret(X), *Reinterpret(Y) ) ) } \
  ElError ElAxpyDist_ ## SIG \
  ( CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY( Axpy( Reinterpret(alpha), *Reinterpret(X), *Reinterpret(Y) ) ) } \
  /* tri(Y) := tri(alpha X + Y) */ \
  ElError ElAxpyTriangle_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY \
    ( AxpyTriangle \
      ( Reinterpret(uplo), Reinterpret(alpha), \
        *Reinterpret(X), *Reinterpret(Y) ) ) } \
  ElError ElAxpyTriangleDist_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY \
    ( AxpyTriangle \
      ( Reinterpret(uplo), Reinterpret(alpha), \
        *Reinterpret(X), *Reinterpret(Y) ) ) } \
  /* B = A */ \
  ElError ElCopy_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *Reinterpret(A), *Reinterpret(B) ) ) } \
  ElError ElCopyDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *Reinterpret(A), *Reinterpret(B) ) ) } \
  /* DiagonalScale */ \
  ElError ElDiagonalScale_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG X ) \
  { EL_TRY( \
      DiagonalScale \
      ( Reinterpret(side), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X) ) ) } \
  ElError ElDiagonalScaleDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( \
      DiagonalScale \
      ( Reinterpret(side), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X) ) ) } \
  /* DiagonalScaleTrapezoid */ \
  ElError ElDiagonalScaleTrapezoid_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG X, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( Reinterpret(side), Reinterpret(uplo), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X), offset ) ) } \
  ElError ElDiagonalScaleTrapezoidDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG X, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( Reinterpret(side), Reinterpret(uplo), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X), offset ) ) } \
  /* Dot product (<A,B>=vec(A)^H vec(B)) */ \
  ElError ElDot_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = Reinterpret(Dot(*Reinterpret(A),*Reinterpret(B))) ) } \
  ElError ElDotDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = Reinterpret(Dot(*Reinterpret(A),*Reinterpret(B))) ) } \
  /* Unconjugated dot product, vec(A)^T vec(B) */ \
  ElError ElDotu_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = Reinterpret(Dotu(*Reinterpret(A),*Reinterpret(B))) ) } \
  ElError ElDotuDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = Reinterpret(Dotu(*Reinterpret(A),*Reinterpret(B))) ) } \
  /* EntrywiseFill */ \
  ElError ElEntrywiseFill_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*fill)() ) \
  { try { \
      auto newFill = [&]() { return Reinterpret(fill()); }; \
      EntrywiseFill( *Reinterpret(A), std::function<T(void)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseFillDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*fill)() ) \
  { try { \
      auto newFill = [&]() { return Reinterpret(fill()); }; \
      EntrywiseFill( *Reinterpret(A), std::function<T(void)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* EntrywiseMap */ \
  ElError ElEntrywiseMap_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return Reinterpret(func(Reinterpret(alpha))); }; \
      EntrywiseMap( *Reinterpret(A), std::function<T(T)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseMapDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return Reinterpret(func(Reinterpret(alpha))); }; \
      EntrywiseMap( *Reinterpret(A), std::function<T(T)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Fill */ \
  ElError ElFill_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *Reinterpret(A), Reinterpret(alpha) ) ) } \
  ElError ElFillDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *Reinterpret(A), Reinterpret(alpha) ) ) } \
  /* Hadamard */ \
  ElError ElHadamard_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG C ) \
  { EL_TRY( Hadamard(*Reinterpret(A),*Reinterpret(B),*Reinterpret(C)) ) } \
  ElError ElHadamardDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG C ) \
  { EL_TRY( Hadamard(*Reinterpret(A),*Reinterpret(B),*Reinterpret(C)) ) } \
  /* Hilbert-Schmidt inner product (same as Dot) */ \
  ElError ElHilbertSchmidt_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = \
      Reinterpret(HilbertSchmidt(*Reinterpret(A),*Reinterpret(B))) ) } \
  ElError ElHilbertSchmidtDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = \
      Reinterpret(HilbertSchmidt(*Reinterpret(A),*Reinterpret(B))) ) } \
  /* IndexDependentFill */ \
  ElError ElIndexDependentFill_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*fill)(ElInt,ElInt) ) \
  { try { \
      auto newFill = [&]( Int i, Int j ) { return Reinterpret(fill(i,j)); }; \
      IndexDependentFill \
      ( *Reinterpret(A), std::function<T(Int,Int)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElIndexDependentFillDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*fill)(ElInt,ElInt) ) \
  { try { \
      auto newFill = [&]( Int i, Int j ) { return Reinterpret(fill(i,j)); }; \
      IndexDependentFill \
      ( *Reinterpret(A), std::function<T(Int,Int)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* IndexDependentMap */ \
  ElError ElIndexDependentMap_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*func)(ElInt,ElInt,CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( Int i, Int j, T alpha ) \
        { return Reinterpret(func(i,j,Reinterpret(alpha))); }; \
      IndexDependentMap \
      ( *Reinterpret(A), std::function<T(Int,Int,T)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElIndexDependentMapDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*func)(ElInt,ElInt,CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( Int i, Int j, T alpha ) \
        { return Reinterpret(func(i,j,Reinterpret(alpha))); }; \
      IndexDependentMap \
      ( *Reinterpret(A), std::function<T(Int,Int,T)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* MakeHermitian */ \
  ElError ElMakeHermitian_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( MakeHermitian( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElMakeHermitianDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( MakeHermitian( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  /* MakeSymmetric */ \
  ElError ElMakeSymmetric_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( MakeSymmetric( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElMakeSymmetricDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( MakeSymmetric( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  /* MakeTrapezoidal */ \
  ElError ElMakeTrapezoidal_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( MakeTrapezoidal( Reinterpret(uplo), *Reinterpret(A), offset ) ) } \
  ElError ElMakeTrapezoidalDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( MakeTrapezoidal( Reinterpret(uplo), *Reinterpret(A), offset ) ) } \
  /* MakeTriangular */ \
  ElError ElMakeTriangular_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( MakeTriangular( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElMakeTriangularDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( MakeTriangular( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  /* TODO: MaxAbs */ \
  /* TODO: MinAbs */ \
  /* TODO: QuasiDiagonalScale */ \
  /* Scale */ \
  ElError ElScale_ ## SIG \
  ( CREFLECT(T) alpha, ElMatrix_ ## SIG A ) \
  { EL_TRY( Scale( Reinterpret(alpha), *Reinterpret(A) ) ) } \
  ElError ElScaleDist_ ## SIG \
  ( CREFLECT(T) alpha, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Scale( Reinterpret(alpha), *Reinterpret(A) ) ) } \
  /* ScaleTrapezoid */ \
  ElError ElScaleTrapezoid_ ## SIG \
  ( CREFLECT(T) alpha, ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( \
      ScaleTrapezoid \
      ( Reinterpret(alpha), Reinterpret(uplo), *Reinterpret(A), offset ) ) } \
  ElError ElScaleTrapezoidDist_ ## SIG \
  ( CREFLECT(T) alpha, ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, \
    ElInt offset ) \
  { EL_TRY( \
      ScaleTrapezoid \
      ( Reinterpret(alpha), Reinterpret(uplo), *Reinterpret(A), offset ) ) } \
  /* SetDiagonal */ \
  ElError ElSetDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( SetDiagonal( *Reinterpret(A), Reinterpret(alpha), offset ) ) } \
  ElError ElSetDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( SetDiagonal( *Reinterpret(A), Reinterpret(alpha), offset ) ) } \
  /* Swap */ \
  ElError ElSwap_ ## SIG \
  ( ElOrientation orientation, ElMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY( \
      Swap( Reinterpret(orientation), *Reinterpret(X), *Reinterpret(Y) ) ) } \
  ElError ElSwapDist_ ## SIG \
  ( ElOrientation orientation, \
    ElDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY( \
      Swap( Reinterpret(orientation), *Reinterpret(X), *Reinterpret(Y) ) ) } \
  ElError ElRowSwap_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( RowSwap( *Reinterpret(A), to, from ) ) } \
  ElError ElRowSwapDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( RowSwap( *Reinterpret(A), to, from ) ) } \
  ElError ElColSwap_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( ColSwap( *Reinterpret(A), to, from ) ) } \
  ElError ElColSwapDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( ColSwap( *Reinterpret(A), to, from ) ) } \
  ElError ElSymmetricSwap_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( SymmetricSwap( Reinterpret(uplo), *Reinterpret(A), to, from ) ) } \
  ElError ElSymmetricSwapDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( SymmetricSwap( Reinterpret(uplo), *Reinterpret(A), to, from ) ) } \
  ElError ElHermitianSwap_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( HermitianSwap( Reinterpret(uplo), *Reinterpret(A), to, from ) ) } \
  ElError ElHermitianSwapDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( HermitianSwap( Reinterpret(uplo), *Reinterpret(A), to, from ) ) } \
  /* TODO: Symmetric2x2Scale */ \
  /* B = A^T */ \
  ElError ElTranspose_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*Reinterpret(A),*Reinterpret(B),false) ) } \
  ElError ElTransposeDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*Reinterpret(A),*Reinterpret(B),false) ) } \
  /* UpdateDiagonal */ \
  ElError ElUpdateDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( UpdateDiagonal( *Reinterpret(A), Reinterpret(alpha), offset ) ) } \
  ElError ElUpdateDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( UpdateDiagonal( *Reinterpret(A), Reinterpret(alpha), offset ) ) } \
  /* Zero */ \
  ElError ElZero_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Zero( *Reinterpret(A) ) ) } \
  ElError ElZeroDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Zero( *Reinterpret(A) ) ) }

#define C_PROTO_NOCOMPLEX(SIG,T) \
  /* Max */ \
  ElError ElMax_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElValueIntPair_ ## SIG *entry ) \
  { EL_TRY( *entry = Reinterpret(Max(*Reinterpret(A))) ) } \
  ElError ElMaxDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElValueIntPair_ ## SIG *entry ) \
  { EL_TRY( *entry = Reinterpret(Max(*Reinterpret(A))) ) } \
  ElError ElSymmetricMax_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    ElValueIntPair_ ## SIG *entry ) \
  { EL_TRY( *entry = \
      Reinterpret(SymmetricMax(Reinterpret(uplo),*Reinterpret(A))) ) } \
  ElError ElSymmetricMaxDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    ElValueIntPair_ ## SIG *entry ) \
  { EL_TRY( *entry = \
      Reinterpret(SymmetricMax(Reinterpret(uplo),*Reinterpret(A))) ) } \
  ElError ElVectorMax_ ## SIG \
  ( ElConstMatrix_ ## SIG x, ElValueInt_ ## SIG *entry ) \
  { EL_TRY( *entry = Reinterpret(VectorMax(*Reinterpret(x))) ) } \
  ElError ElVectorMaxDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG x, ElValueInt_ ## SIG *entry ) \
  { EL_TRY( *entry = Reinterpret(VectorMax(*Reinterpret(x))) ) } \
  /* TODO: DiagonalMax */ \
  /* Min */ \
  ElError ElMin_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElValueIntPair_ ## SIG *entry ) \
  { EL_TRY( *entry = Reinterpret(Min(*Reinterpret(A))) ) } \
  ElError ElMinDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElValueIntPair_ ## SIG *entry ) \
  { EL_TRY( *entry = Reinterpret(Min(*Reinterpret(A))) ) } \
  ElError ElSymmetricMin_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    ElValueIntPair_ ## SIG *entry ) \
  { EL_TRY( *entry = \
      Reinterpret(SymmetricMin(Reinterpret(uplo),*Reinterpret(A))) ) } \
  ElError ElSymmetricMinDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    ElValueIntPair_ ## SIG *entry ) \
  { EL_TRY( *entry = \
      Reinterpret(SymmetricMin(Reinterpret(uplo),*Reinterpret(A))) ) } \
  ElError ElVectorMin_ ## SIG \
  ( ElConstMatrix_ ## SIG x, ElValueInt_ ## SIG *entry ) \
  { EL_TRY( *entry = Reinterpret(VectorMin(*Reinterpret(x))) ) } \
  ElError ElVectorMinDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG x, ElValueInt_ ## SIG *entry ) \
  { EL_TRY( *entry = Reinterpret(VectorMin(*Reinterpret(x))) ) } \
  /* TODO: DiagonalMin */

#define C_PROTO_NOINT(SIG,T) \
  /* DiagonalSolve */ \
  ElError ElDiagonalSolve_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG X ) \
  { EL_TRY( \
      DiagonalSolve \
      ( Reinterpret(side), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X) ) ) } \
  ElError ElDiagonalSolveDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( \
      DiagonalSolve \
      ( Reinterpret(side), Reinterpret(orientation), \
        *Reinterpret(d), *Reinterpret(X) ) ) } \
  /* Nrm2 (same as FrobeniusNorm) */ \
  ElError ElNrm2_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<T> *gamma ) \
  { EL_TRY( *gamma = Nrm2(*Reinterpret(A)) ) } \
  ElError ElNrm2Dist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<T> *gamma ) \
  { EL_TRY( *gamma = Nrm2(*Reinterpret(A)) ) } \
  /* TODO: QuasiDiagonalSolve */ \
  /* TODO: Symmetric2x2Inv */ \
  /* TODO: Symmetric2x2Solve */

#define C_PROTO_INT(SIG,T) \
  C_PROTO_BASE(SIG,T) \
  C_PROTO_NOCOMPLEX(SIG,T)

#define C_PROTO_REAL(SIG,T) \
  C_PROTO_BASE(SIG,T) \
  C_PROTO_NOCOMPLEX(SIG,T) \
  C_PROTO_NOINT(SIG,T)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO_BASE(SIG,T) \
  C_PROTO_NOINT(SIG,T) \
  /* Conjugate */ \
  ElError ElConjugate_ ## SIG( ElMatrix_ ## SIG A ) \
  { EL_TRY( Conjugate( *Reinterpret(A) ) ) } \
  ElError ElConjugateDist_ ## SIG( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Conjugate( *Reinterpret(A) ) ) } \
  /* MakeReal */ \
  ElError ElMakeReal_ ## SIG \
  ( ElMatrix_ ## SIG A ) { EL_TRY( MakeReal( *Reinterpret(A) ) ) } \
  ElError ElMakeRealDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A ) { EL_TRY( MakeReal( *Reinterpret(A) ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
