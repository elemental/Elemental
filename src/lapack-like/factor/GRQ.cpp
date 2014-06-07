/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F> 
void GRQ( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    Matrix<F> tA;
    Matrix<Base<F>> dA;
    RQ( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    QR( B );
}

template<typename F> 
void GRQ( DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> tA(g);
    DistMatrix<Base<F>,MD,STAR> dA(g);
    RQ( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    QR( B );
}

template<typename F> 
void GRQ
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, 
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    RQ( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    QR( B, tB, dB );
}

template<typename F> 
void GRQ
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& tA, DistMatrix<Base<F>,MD,STAR>& dA,
  DistMatrix<F>& B, DistMatrix<F,MD,STAR>& tB, DistMatrix<Base<F>,MD,STAR>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    RQ( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    QR( B, tB, dB );
}

#define PROTO(F) \
  template void GRQ( Matrix<F>& A, Matrix<F>& B ); \
  template void GRQ( DistMatrix<F>& A, DistMatrix<F>& B ); \
  template void GRQ \
  ( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, \
    Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB ); \
  template void GRQ \
  ( DistMatrix<F>& A, \
    DistMatrix<F,MD,STAR>& tA, DistMatrix<Base<F>,MD,STAR>& dA, \
    DistMatrix<F>& B, \
    DistMatrix<F,MD,STAR>& tB, DistMatrix<Base<F>,MD,STAR>& dB );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
