/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LQ_EXPLICIT_HPP
#define ELEM_LAPACK_LQ_EXPLICIT_HPP

#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/lapack-like/LQ.hpp"
#include "elemental/matrices/Identity.hpp"

namespace elem {
namespace lq {

template<typename F>
inline void
Explicit( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Explicit"))
    Matrix<F> t;
    LQ( A, t );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<F> Q;
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, t, Q );
    A = Q;
}

template<typename F>
inline void
Explicit( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t( g );
    LQ( A, t );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<F> Q( g );
    Q.AlignWith( A );
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, t, Q );
    A = Q;
}

template<typename F>
inline void
Explicit( Matrix<F>& L, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Explicit"))
    Matrix<F> t;
    LQ( A, t );
    L = A;
    MakeTriangular( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<F> Q;
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, t, Q );
    A = Q;
}

template<typename F>
inline void
Explicit( DistMatrix<F>& L, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t( g );
    LQ( A, t );
    L = A;
    MakeTriangular( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<F> Q( g );
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, t, Q );
    A = Q;
}

} // namespace lq
} // namespace elem

#endif // ifndef ELEM_LAPACK_LQ_EXPLICIT_HPP
