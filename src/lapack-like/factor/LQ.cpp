/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./LQ/ApplyQ.hpp"
#include "./LQ/Householder.hpp"
#include "./LQ/SolveAfter.hpp"
#include "./LQ/Explicit.hpp"

namespace El {

template<typename F> 
void LQ( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LQ"))
    lq::Householder( A );
}

template<typename F> 
void LQ( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LQ"))
    lq::Householder( A );
}

template<typename F> 
void LQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("LQ"))
    lq::Householder( A, t, d );
}

template<typename F> 
void LQ
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(CallStackEntry cse("LQ"))
    lq::Householder( A, t, d );
}

// Variants which perform (Businger-Golub) row-pivoting
// ====================================================
// TODO

#define PROTO(F) \
  template void LQ( Matrix<F>& A ); \
  template void LQ( DistMatrix<F>& A ); \
  template void LQ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d ); \
  template void LQ \
  ( DistMatrix<F>& A, \
    DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d ); \
  template void lq::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, Matrix<F>& B ); \
  template void lq::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, \
    const DistMatrix<Base<F>,MD,STAR>& d, DistMatrix<F>& B ); \
  template void lq::Explicit( Matrix<F>& A ); \
  template void lq::Explicit( DistMatrix<F>& A ); \
  template void lq::Explicit( Matrix<F>& L, Matrix<F>& A ); \
  template void lq::Explicit( DistMatrix<F>& L, DistMatrix<F>& A ); \
  template void lq::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void lq::SolveAfter \
  ( Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, \
    const DistMatrix<Base<F>,MD,STAR>& d, const DistMatrix<F>& B, \
          DistMatrix<F>& X );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
