/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and 
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./LDL/dense/Var3.hpp"

#include "./LDL/dense/Pivoted.hpp"

#include "./LDL/dense/MultiplyAfter.hpp"
#include "./LDL/dense/SolveAfter.hpp"

#include "./LDL/dense/Inertia.hpp"

#include "./LDL/sparse/numeric/Process.hpp"

namespace El {

// Dense
// =====

// Unpivoted
// ---------
template<typename F>
void LDL( Matrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CSE cse("LDL"))
    ldl::Var3( A, conjugate );
}

template<typename F>
void LDL( ElementalMatrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CSE cse("LDL"))
    ldl::Var3( A, conjugate );
}

template<typename F>
void LDL( DistMatrix<F,STAR,STAR>& A, bool conjugate )
{ LDL( A.Matrix(), conjugate ); }

// Pivoted
// -------
template<typename F>
void LDL
( Matrix<F>& A,
  Matrix<F>& dSub, 
  Permutation& p,
  bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LDL"))
    ldl::Pivoted( A, dSub, p, conjugate, ctrl );
}

template<typename F>
void LDL
( ElementalMatrix<F>& A,
  ElementalMatrix<F>& dSub, 
  DistPermutation& p,
  bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl ) 
{
    DEBUG_ONLY(CSE cse("LDL"))
    ldl::Pivoted( A, dSub, p, conjugate, ctrl );
}

// Sparse
// ======
template<typename F>
void LDL
( const ldl::NodeInfo& info,
        ldl::Front<F>& front,
  LDLFrontType newType )
{
    DEBUG_ONLY(CSE cse("LDL"))
    if( !Unfactored(front.type) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( front, SYMM_2D );

    // Perform the initial factorization
    ldl::Process( info, front, InitialFactorType(newType) );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( front, newType );
}

template<typename F>
void LDL
( const ldl::DistNodeInfo& info,
        ldl::DistFront<F>& front, 
  LDLFrontType newType )
{
    DEBUG_ONLY(CSE cse("LDL"))
    if( !Unfactored(front.type) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( front, SYMM_2D );

    // Perform the initial factorization
    ldl::Process( info, front, InitialFactorType(newType) );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( front, newType );
}

#define PROTO_BASE(F) \
  template void LDL( Matrix<F>& A, bool conjugate ); \
  template void LDL( ElementalMatrix<F>& A, bool conjugate ); \
  template void LDL( DistMatrix<F,STAR,STAR>& A, bool conjugate ); \
  template void LDL \
  ( Matrix<F>& A, \
    Matrix<F>& dSub, \
    Permutation& p, \
    bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void LDL \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& dSub, \
    DistPermutation& p, \
    bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template InertiaType ldl::Inertia \
  ( const Matrix<Base<F>>& d, \
    const Matrix<F>& dSub ); \
  template InertiaType ldl::Inertia \
  ( const ElementalMatrix<Base<F>>& d, \
    const ElementalMatrix<F>& dSub ); \
  template void ldl::MultiplyAfter \
  ( const Matrix<F>& A, \
          Matrix<F>& B, \
    bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const ElementalMatrix<F>& A, \
          ElementalMatrix<F>& B, \
    bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const Matrix<F>& A, \
    const Matrix<F>& dSub, \
    const Permutation& p, \
          Matrix<F>& B, \
    bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& dSub, \
    const DistPermutation& p, \
          ElementalMatrix<F>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const Matrix<F>& A, \
          Matrix<F>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const ElementalMatrix<F>& A, \
          ElementalMatrix<F>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const Matrix<F>& A, \
    const Matrix<F>& dSub, \
    const Permutation& p, \
          Matrix<F>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& dSub, \
    const DistPermutation& p, \
          ElementalMatrix<F>& B, \
     bool conjugated );

#define PROTO(F) \
  PROTO_BASE(F) \
  template void LDL \
  ( const ldl::NodeInfo& info, \
          ldl::Front<F>& front, \
    LDLFrontType newType ); \
  template void LDL \
  ( const ldl::DistNodeInfo& info, \
          ldl::DistFront<F>& front, \
    LDLFrontType newType );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
