/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LQ_EXPLICIT_HPP
#define ELEM_LQ_EXPLICIT_HPP

#include ELEM_MAKETRIANGULAR_INC
#include ELEM_LQ_INC
#include ELEM_IDENTITY_INC

namespace elem {
namespace lq {

template<typename F>
inline void
Explicit( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    LQ( A, t, d );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<F> Q;
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, t, d, Q );
    A = Q;
}

template<typename F>
inline void
Explicit( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    LQ( A, t, d );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<F> Q(g);
    Q.AlignWith( A );
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, t, d, Q );
    A = Q;
}

template<typename F>
inline void
Explicit( Matrix<F>& L, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Explicit"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    LQ( A, t, d );
    Matrix<F> AL, AR;
    PartitionRight( A, AL, AR, Min(A.Height(),A.Width()) );
    L = AL;
    MakeTriangular( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<F> Q;
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, t, d, Q );
    A = Q;
}

template<typename F>
inline void
Explicit( DistMatrix<F>& L, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::Explicit"))
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    LQ( A, t, d );
    DistMatrix<F> AL(g), AR(g);
    PartitionRight( A, AL, AR, Min(A.Height(),A.Width()) );
    L = AL;
    MakeTriangular( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<F> Q(g);
    Identity( Q, A.Height(), A.Width() );
    lq::ApplyQ( RIGHT, NORMAL, A, t, d, Q );
    A = Q;
}

} // namespace lq
} // namespace elem

#endif // ifndef ELEM_LQ_EXPLICIT_HPP
