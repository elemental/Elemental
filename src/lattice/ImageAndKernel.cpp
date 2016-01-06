/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void LatticeImageAndKernel
( Matrix<F>& B,
  Matrix<F>& M,
  Matrix<F>& K,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LatticeImageAndKernel"))

    // NOTE: UInv and R don't actually need to be formed...but deciding on 
    //       the best interface is somewhat tricky
    Matrix<F> U, UInv, R;
    auto info = LLL( B, U, UInv, R, ctrl );
    const Int rank = info.rank;
    const Int n = B.Width();
    M = B(ALL,IR(0,rank));
    K = U(ALL,IR(rank,n));

    // Reduce the columns of U that corresponded to the kernel
    LLL( K, ctrl );

    // Rather than explicitly inverting the Gram matrix of the kernel basis
    // as suggested by Cohen, we can simply solve a least squares problem
    //
    // NOTE: 'R' is reused for the least squares solution
    LeastSquares( NORMAL, K, M, R );
    Round( R );
    Gemm( NORMAL, NORMAL, F(-1), K, R, F(1), M );
}

template<typename F>
void LatticeKernel
( Matrix<F>& B,
  Matrix<F>& K,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LatticeKernel"))
    // TODO: Take the shortcuts suggested in Algorithm 2.7.2 of Henri Cohen's
    //       "A course in computational algebraic number theory".
    Matrix<F> M;
    LatticeImageAndKernel( B, M, K, ctrl );
}

#define PROTO(F) \
  template void LatticeImageAndKernel \
  ( Matrix<F>& B, \
    Matrix<F>& M, \
    Matrix<F>& K, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template void LatticeKernel \
  ( Matrix<F>& B, \
    Matrix<F>& K, \
    const LLLCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
