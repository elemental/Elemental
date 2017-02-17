/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2016, Ron Estrin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
void LatticeImageAndKernel
( const Matrix<Field>& B,
        Matrix<Field>& M,
        Matrix<Field>& K,
  const LLLCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE

    Matrix<Field> BCopy( B );
    Matrix<Field> U, R;
    auto info = LLL( BCopy, U, R, ctrl );
    const Int rank = info.rank;
    const Int n = B.Width();
    M = BCopy(ALL,IR(0,rank));
    K = U(ALL,IR(rank,n));

    // Reduce the columns of U that corresponded to the kernel
    LLL( K, ctrl );

    // Rather than explicitly inverting the Gram matrix of the kernel basis
    // as suggested by Cohen, we can simply solve a least squares problem
    //
    // NOTE: 'R' is reused for the least squares solution
    // TODO(poulson): Support other options than just "Babai rounding", e.g.,
    // Nulling and Cancelling (with optimal ordering)
    LeastSquares( NORMAL, K, M, R );
    Round( R );
    Gemm( NORMAL, NORMAL, Field(-1), K, R, Field(1), M );
}

template<typename Field>
void LatticeImage
( const Matrix<Field>& B,
        Matrix<Field>& M,
  const LLLCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    // TODO(poulson): Avoid computing the kernel?
    Matrix<Field> K;
    LatticeImageAndKernel( B, M, K, ctrl );
}

template<typename Field>
void LatticeKernel
( const Matrix<Field>& B,
        Matrix<Field>& K,
  const LLLCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    // TODO(poulson): Take the shortcuts suggested in Algorithm 2.7.2 of
    // Henri Cohen's
    //       "A course in computational algebraic number theory".

    Matrix<Field> BCopy( B );
    Matrix<Field> U, R;
    auto info = LLL( BCopy, U, R, ctrl );
    const Int rank = info.rank;
    const Int n = B.Width();
    K = U(ALL,IR(rank,n));

    // Reduce the columns of U that corresponded to the kernel
    LLL( K, ctrl );
}

#define PROTO(Field) \
  template void LatticeImageAndKernel \
  ( const Matrix<Field>& B, \
          Matrix<Field>& M, \
          Matrix<Field>& K, \
    const LLLCtrl<Base<Field>>& ctrl ); \
  template void LatticeImage \
  ( const Matrix<Field>& B, \
          Matrix<Field>& M, \
    const LLLCtrl<Base<Field>>& ctrl ); \
  template void LatticeKernel \
  ( const Matrix<Field>& B, \
          Matrix<Field>& K, \
    const LLLCtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
