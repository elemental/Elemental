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

#define C_PROTO(SIG,T) \
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
  /* Conjugate */ \
  ElError ElConjugateMatrix_ ## SIG( ElMatrix_ ## SIG A ) \
  { EL_TRY( Conjugate( *Reinterpret(A) ) ) } \
  ElError ElConjugateDistMatrix_ ## SIG( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Conjugate( *Reinterpret(A) ) ) } \
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
  /* TODO: DiagonalScaleTrapezoid */ \
  /* B = A^T */ \
  ElError ElTransposeMatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*Reinterpret(A),*Reinterpret(B),false) ) } \
  ElError ElTransposeDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*Reinterpret(A),*Reinterpret(B),false) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
