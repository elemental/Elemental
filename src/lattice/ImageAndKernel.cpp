/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void LatticeImageAndKernel
( const Matrix<F>& B,
        Matrix<F>& M,
        Matrix<F>& K,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LatticeImageAndKernel"))

    Matrix<F> BCopy( B );
    Matrix<F> U, R;
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
    // TODO: Support other options than just "Babai rounding", e.g.,
    //       Nulling and Cancelling (with optimal ordering)
    LeastSquares( NORMAL, K, M, R );
    Round( R );
    Gemm( NORMAL, NORMAL, F(-1), K, R, F(1), M );
}

template<typename F>
void LatticeImage
( const Matrix<F>& B,
        Matrix<F>& M,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LatticeImage"))
    // TODO: Avoid computing the kernel?
    Matrix<F> K;
    LatticeImageAndKernel( B, M, K, ctrl );
}

template<typename F>
void LatticeKernel
( const Matrix<F>& B,
        Matrix<F>& K,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LatticeKernel"))
    // TODO: Take the shortcuts suggested in Algorithm 2.7.2 of Henri Cohen's
    //       "A course in computational algebraic number theory".

    Matrix<F> BCopy( B );
    Matrix<F> U, R;
    auto info = LLL( BCopy, U, R, ctrl );
    const Int rank = info.rank;
    const Int n = B.Width();
    K = U(ALL,IR(rank,n));

    // Reduce the columns of U that corresponded to the kernel
    LLL( K, ctrl );
}

#define PROTO(F) \
  template void LatticeImageAndKernel \
  ( const Matrix<F>& B, \
          Matrix<F>& M, \
          Matrix<F>& K, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template void LatticeImage \
  ( const Matrix<F>& B, \
          Matrix<F>& M, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template void LatticeKernel \
  ( const Matrix<F>& B, \
          Matrix<F>& K, \
    const LLLCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
