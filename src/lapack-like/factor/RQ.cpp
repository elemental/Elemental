/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./RQ/ApplyQ.hpp"
#include "./RQ/Cholesky.hpp"
#include "./RQ/Householder.hpp"
#include "./RQ/SolveAfter.hpp"

namespace El {

template<typename F> 
void RQ( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("RQ"))
    rq::Householder( A );
}

template<typename F> 
void RQ( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("RQ"))
    rq::Householder( A );
}

template<typename F> 
void RQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("RQ"))
    rq::Householder( A, t, d );
}

template<typename F> 
void RQ
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(CallStackEntry cse("RQ"))
    rq::Householder( A, t, d );
}

// Variants which perform (Businger-Golub) row-pivoting
// ====================================================
// TODO

#define PROTO(F) \
  template void RQ( Matrix<F>& A ); \
  template void RQ( DistMatrix<F>& A ); \
  template void RQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d ); \
  template void RQ \
  ( DistMatrix<F>& A, \
    DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d ); \
  template void rq::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, Matrix<F>& B ); \
  template void rq::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, \
    const DistMatrix<Base<F>,MD,STAR>& d, DistMatrix<F>& B ); \
  template void rq::Cholesky( Matrix<F>& A, Matrix<F>& R ); \
  template void rq::Cholesky \
  ( DistMatrix<F,STAR,VR>& A, DistMatrix<F,STAR,STAR>& R ); \
  template void rq::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void rq::SolveAfter \
  ( Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, \
    const DistMatrix<Base<F>,MD,STAR>& d, const DistMatrix<F>& B, \
          DistMatrix<F>& X );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
