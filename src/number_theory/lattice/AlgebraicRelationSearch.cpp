/*
   Copyright (c) 2009-2016, Jack Poulson, 2016, Ron Estrin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// TODO: A version which pushes for small coefficients for the large degrees
//       by modifying the identity matrix in the upper portion of B

template<typename F>
Int AlgebraicRelationSearch
( F alpha,
  Int n,
  Base<F> NSqrt,
  Matrix<F>& B,
  Matrix<F>& U,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = n+1;

    Identity( B, m, n );
    auto bLastRow = B( IR(m-1), ALL );
    for( Int j=0; j<n; ++j )
        bLastRow(0,j) = Pow(alpha,Real(j));
    Scale( NSqrt, bLastRow );

    Matrix<F> R;
    auto info = LLL( B, U, R, ctrl );

    return info.nullity;
}

#define PROTO(F) \
  template Int AlgebraicRelationSearch \
  ( F alpha, \
    Int n, \
    Base<F> NSqrt, \
    Matrix<F>& B, \
    Matrix<F>& U, \
    const LLLCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
