/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GRQ_HPP
#define ELEM_GRQ_HPP

#include ELEM_QR_INC
#include ELEM_RQ_INC

namespace elem {

template<typename F> 
inline void
GRQ( Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    Matrix<F> tA;
    Matrix<Base<F>> dA;
    rq::Householder( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    qr::Householder( B );
}

template<typename F> 
inline void
GRQ( DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    const Grid& g = A.Grid();
    DistMatrix<F> tA(g);
    DistMatrix<Base<F>> dA(g);
    rq::Householder( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    MakeTrapezoidal( UPPER, A, Min(A.Height(),A.Width()) );
    qr::Householder( B );
}

template<typename F> 
inline void
GRQ
( Matrix<F>& A, Matrix<F>& tA, Matrix<Base<F>>& dA, 
  Matrix<F>& B, Matrix<F>& tB, Matrix<Base<F>>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    rq::Householder( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    qr::Householder( B, tB, dB );
}

template<typename F> 
inline void
GRQ
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& tA, DistMatrix<Base<F>,MD,STAR>& dA,
  DistMatrix<F>& B, DistMatrix<F,MD,STAR>& tB, DistMatrix<Base<F>,MD,STAR>& dB )
{
    DEBUG_ONLY(CallStackEntry cse("GRQ"))
    rq::Householder( A, tA, dA );
    rq::ApplyQ( RIGHT, ADJOINT, A, tA, dA, B );
    qr::Householder( B, tB, dB );
}

} // namespace elem

#endif // ifndef ELEM_GRQ_HPP
