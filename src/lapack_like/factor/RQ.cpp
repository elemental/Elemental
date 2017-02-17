/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./RQ/ApplyQ.hpp"
#include "./RQ/Householder.hpp"

#include "./RQ/SolveAfter.hpp"

#include "./RQ/Explicit.hpp"

#include "./RQ/Cholesky.hpp"

namespace El {

template<typename F>
void RQ
( Matrix<F>& A,
  Matrix<F>& householderScalars,
  Matrix<Base<F>>& signature )
{
    EL_DEBUG_CSE
    rq::Householder( A, householderScalars, signature );
}

template<typename F>
void RQ
( AbstractDistMatrix<F>& A,
  AbstractDistMatrix<F>& householderScalars,
  AbstractDistMatrix<Base<F>>& signature )
{
    EL_DEBUG_CSE
    rq::Householder( A, householderScalars, signature );
}

// Variants which perform (Businger-Golub) row-pivoting
// ====================================================
// TODO

#define PROTO(F) \
  template void RQ \
  ( Matrix<F>& A, \
    Matrix<F>& householderScalars, \
    Matrix<Base<F>>& signature ); \
  template void RQ \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& householderScalars, \
    AbstractDistMatrix<Base<F>>& signature ); \
  template void rq::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
    const Matrix<Base<F>>& signature, \
          Matrix<F>& B ); \
  template void rq::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
    const AbstractDistMatrix<Base<F>>& signature, \
          AbstractDistMatrix<F>& B ); \
  template void rq::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
    const Matrix<Base<F>>& signature, \
    const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void rq::SolveAfter \
  ( Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
    const AbstractDistMatrix<Base<F>>& signature, \
    const AbstractDistMatrix<F>& B, \
          AbstractDistMatrix<F>& X ); \
  template void rq::Cholesky( Matrix<F>& A, Matrix<F>& R ); \
  template void rq::Cholesky \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R ); \
  template void rq::ExplicitTriang( Matrix<F>& A ); \
  template void rq::ExplicitTriang( AbstractDistMatrix<F>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
