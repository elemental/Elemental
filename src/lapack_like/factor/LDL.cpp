/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013 Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2014 Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./LDL/dense/Var3.hpp"

#include "./LDL/dense/Pivoted.hpp"

#include "./LDL/dense/MultiplyAfter.hpp"
#include "./LDL/dense/SolveAfter.hpp"

#include "./LDL/dense/Inertia.hpp"

namespace El {

// Dense
// =====

// Unpivoted
// ---------
template<typename Field>
void LDL( Matrix<Field>& A, bool conjugate )
{
    EL_DEBUG_CSE
    ldl::Var3( A, conjugate );
}

template<typename Field>
void LDL( AbstractDistMatrix<Field>& A, bool conjugate )
{
    EL_DEBUG_CSE
    ldl::Var3( A, conjugate );
}

template<typename Field>
void LDL( DistMatrix<Field,STAR,STAR>& A, bool conjugate )
{ LDL( A.Matrix(), conjugate ); }

// Pivoted
// -------
template<typename Field>
void LDL
( Matrix<Field>& A,
  Matrix<Field>& dSub,
  Permutation& p,
  bool conjugate,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    ldl::Pivoted( A, dSub, p, conjugate, ctrl );
}

template<typename Field>
void LDL
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& dSub,
  DistPermutation& p,
  bool conjugate,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    ldl::Pivoted( A, dSub, p, conjugate, ctrl );
}

#define PROTO(Field) \
  template void LDL( Matrix<Field>& A, bool conjugate ); \
  template void LDL( AbstractDistMatrix<Field>& A, bool conjugate ); \
  template void LDL( DistMatrix<Field,STAR,STAR>& A, bool conjugate ); \
  template void LDL \
  ( Matrix<Field>& A, \
    Matrix<Field>& dSub, \
    Permutation& p, \
    bool conjugate, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void LDL \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Field>& dSub, \
    DistPermutation& p, \
    bool conjugate, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template InertiaType ldl::Inertia \
  ( const Matrix<Base<Field>>& d, \
    const Matrix<Field>& dSub ); \
  template InertiaType ldl::Inertia \
  ( const AbstractDistMatrix<Base<Field>>& d, \
    const AbstractDistMatrix<Field>& dSub ); \
  template void ldl::MultiplyAfter \
  ( const Matrix<Field>& A, \
          Matrix<Field>& B, \
    bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& B, \
    bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const Matrix<Field>& A, \
    const Matrix<Field>& dSub, \
    const Permutation& p, \
          Matrix<Field>& B, \
    bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const AbstractDistMatrix<Field>& A, \
    const AbstractDistMatrix<Field>& dSub, \
    const DistPermutation& p, \
          AbstractDistMatrix<Field>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const Matrix<Field>& A, \
          Matrix<Field>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const Matrix<Field>& A, \
    const Matrix<Field>& dSub, \
    const Permutation& p, \
          Matrix<Field>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const AbstractDistMatrix<Field>& A, \
    const AbstractDistMatrix<Field>& dSub, \
    const DistPermutation& p, \
          AbstractDistMatrix<Field>& B, \
     bool conjugated );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
