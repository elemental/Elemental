/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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
void GQR( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& BPre )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))

    const Grid& g = APre.Grid();
    DistMatrix<F> A(g), B(g);
    Copy( APre, A, READ_WRITE_PROXY );
    Copy( BPre, B, READ_WRITE_PROXY );

    DistMatrix<F,MD,STAR> tA(g);
    DistMatrix<Base<F>,MD,STAR> dA(g);
    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    MakeTriangular( UPPER, A );
    RQ( B );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( B, BPre, RESTORE_READ_WRITE_PROXY );
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
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<F>& tA, AbstractDistMatrix<Base<F>>& dA,
  AbstractDistMatrix<F>& BPre, 
  AbstractDistMatrix<F>& tB, AbstractDistMatrix<Base<F>>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))

    const Grid& g = APre.Grid();
    DistMatrix<F> A(g), B(g);
    Copy( APre, A, READ_WRITE_PROXY );
    Copy( BPre, B, READ_WRITE_PROXY );

    QR( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    RQ( B, tB, dB );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( B, BPre, RESTORE_READ_WRITE_PROXY );
}

#define PROTO(F) \
  template void GQR( Matrix<F>& A, Matrix<F>& B ); \
  template void GQR( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B ); \
  template void GQR \
  ( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, \
    Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB ); \
  template void GQR \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& tA, AbstractDistMatrix<Base<F>>& dA, \
    AbstractDistMatrix<F>& B, \
    AbstractDistMatrix<F>& tB, AbstractDistMatrix<Base<F>>& dB );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
