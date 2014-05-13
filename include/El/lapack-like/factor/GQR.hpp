/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_GQR_HPP
#define EL_GQR_HPP

#include EL_QR_INC
#include EL_RQ_INC

namespace El {

template<typename F> 
inline void
GQR( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))
    Matrix<F> tA;
    Matrix<Base<F>> dA;
    qr::Householder( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    MakeTriangular( UPPER, A );
    rq::Householder( B );
}

template<typename F> 
inline void
GQR( DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))
    const Grid& g = A.Grid();
    DistMatrix<F> tA(g);
    DistMatrix<Base<F>> dA(g);
    qr::Householder( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    MakeTriangular( UPPER, A );
    rq::Householder( B );
}

template<typename F> 
inline void
GQR
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, 
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))
    qr::Householder( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    rq::Householder( B, tB, dB );
}

template<typename F> 
inline void
GQR
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& tA, DistMatrix<Base<F>,MD,STAR>& dA,
  DistMatrix<F>& B, DistMatrix<F,MD,STAR>& tB, DistMatrix<Base<F>,MD,STAR>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GQR"))
    qr::Householder( A, tA, dA );
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, B );
    rq::Householder( B, tB, dB );
}

} // namespace El

#endif // ifndef EL_GQR_HPP
