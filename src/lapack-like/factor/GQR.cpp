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
void GQR( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))
    Matrix<F> tA;
    Matrix<Base<F>> dA;
    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    MakeTriangular( UPPER, A );
    RQ( B );
}

template<typename F> 
void GQR( DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> tA(g);
    DistMatrix<Base<F>,MD,STAR> dA(g);
    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    MakeTriangular( UPPER, A );
    RQ( B );
}

template<typename F> 
void GQR
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, 
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))
    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    RQ( B, tB, dB );
}

template<typename F> 
void GQR
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& tA, DistMatrix<Base<F>,MD,STAR>& dA,
  DistMatrix<F>& B, DistMatrix<F,MD,STAR>& tB, DistMatrix<Base<F>,MD,STAR>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))
    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    RQ( B, tB, dB );
}

#define PROTO(F) \
  template void GQR( Matrix<F>& A, Matrix<F>& B ); \
  template void GQR( DistMatrix<F>& A, DistMatrix<F>& B ); \
  template void GQR \
  ( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, \
    Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB ); \
  template void GQR \
  ( DistMatrix<F>& A, \
    DistMatrix<F,MD,STAR>& tA, DistMatrix<Base<F>,MD,STAR>& dA, \
    DistMatrix<F>& B, \
    DistMatrix<F,MD,STAR>& tB, DistMatrix<Base<F>,MD,STAR>& dB );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
