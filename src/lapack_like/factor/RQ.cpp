/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./RQ/ApplyQ.hpp"
#include "./RQ/Householder.hpp"

#include "./RQ/SolveAfter.hpp"

#include "./RQ/Explicit.hpp"

#include "./RQ/Cholesky.hpp"

namespace El {

template<typename F> 
void RQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CSE cse("RQ"))
    rq::Householder( A, t, d );
}

template<typename F> 
void RQ
( ElementalMatrix<F>& A, ElementalMatrix<F>& t, 
  ElementalMatrix<Base<F>>& d )
{
    DEBUG_ONLY(CSE cse("RQ"))
    rq::Householder( A, t, d );
}

// Variants which perform (Businger-Golub) row-pivoting
// ====================================================
// TODO

#define PROTO(F) \
  template void RQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d ); \
  template void RQ \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& t, ElementalMatrix<Base<F>>& d ); \
  template void rq::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, Matrix<F>& B ); \
  template void rq::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<F>& A, const ElementalMatrix<F>& t, \
    const ElementalMatrix<Base<F>>& d, ElementalMatrix<F>& B ); \
  template void rq::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void rq::SolveAfter \
  ( Orientation orientation, \
    const ElementalMatrix<F>& A, const ElementalMatrix<F>& t, \
    const ElementalMatrix<Base<F>>& d, const ElementalMatrix<F>& B, \
          ElementalMatrix<F>& X ); \
  template void rq::Cholesky( Matrix<F>& A, Matrix<F>& R ); \
  template void rq::Cholesky \
  ( ElementalMatrix<F>& A, ElementalMatrix<F>& R ); \
  template void rq::ExplicitTriang( Matrix<F>& A ); \
  template void rq::ExplicitTriang( ElementalMatrix<F>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
