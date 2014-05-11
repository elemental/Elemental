/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QR_HPP
#define ELEM_QR_HPP

#include "./QR/ApplyQ.hpp"
#include "./QR/BusingerGolub.hpp"
#include "./QR/Cholesky.hpp"
#include "./QR/Householder.hpp"
#include "./QR/SolveAfter.hpp"
#include "./QR/Explicit.hpp"
#include "./QR/TS.hpp"

namespace elem {

template<typename F> 
inline void
QR( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A );
}

template<typename F> 
inline void
QR( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A );
}

template<typename F> 
inline void
QR( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A, t, d );
}

template<typename F> 
inline void
QR( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A, t, d );
}

// Variants which perform (Businger-Golub) column-pivoting
// =======================================================

template<typename F> 
inline void
QR( Matrix<F>& A, Matrix<Int>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::BusingerGolub( A, pPerm );
}

template<typename F,Dist UPerm> 
inline void
QR( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::BusingerGolub( A, pPerm );
}

template<typename F> 
inline void
QR( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d, Matrix<Int>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::BusingerGolub( A, t, d, pPerm );
}

template<typename F,Dist UPerm> 
inline void
QR( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, 
    DistMatrix<Base<F>,MD,STAR>& d, DistMatrix<Int,UPerm,STAR>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::BusingerGolub( A, t, d, pPerm );
}

} // namespace elem

#endif // ifndef ELEM_QR_HPP
