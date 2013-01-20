/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_EXPLICITLQ_HPP
#define LAPACK_EXPLICITLQ_HPP

#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/LQ.hpp"

namespace elem {

namespace internal {

template<typename Real>
inline void
ExplicitLQHelper( Matrix<Real>& A )
{
    LQ( A );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Real> Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( DistMatrix<Real>& A )
{
    LQ( A );

    // TODO: Replace this with an in-place expansion of Q
    const Grid& g = A.Grid();
    DistMatrix<Real> Q( g );
    Q.AlignWith( A );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( Matrix<Complex<Real> >& A )
{
    Matrix<Complex<Real> > t;
    LQ( A, t );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Complex<Real> > Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( DistMatrix<Complex<Real> >& A )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    LQ( A, t );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<Complex<Real> > Q( g );
    Q.AlignWith( A );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( Matrix<Real>& L, Matrix<Real>& A )
{
    LQ( A );
    L = A;
    MakeTrapezoidal( LEFT, LOWER, 0, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Real> Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( DistMatrix<Real>& L, DistMatrix<Real>& A )
{
    LQ( A );
    L = A;
    MakeTrapezoidal( LEFT, LOWER, 0, L );

    // TODO: Replace this with an in-place expansion of Q
    const Grid& g = A.Grid();
    DistMatrix<Real> Q( g );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper( Matrix<Complex<Real> >& L, Matrix<Complex<Real> >& A )
{
    Matrix<Complex<Real> > t;
    LQ( A, t );
    L = A;
    MakeTrapezoidal( LEFT, LOWER, 0, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Complex<Real> > Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitLQHelper
( DistMatrix<Complex<Real> >& L, DistMatrix<Complex<Real> >& A )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    LQ( A, t );
    L = A;
    MakeTrapezoidal( LEFT, LOWER, 0, L );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<Complex<Real> > Q( g );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

} // namespace internal

template<typename F> 
inline void
ExplicitLQ( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitLQ");
#endif
    internal::ExplicitLQHelper( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitLQ( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitLQ");
#endif
    internal::ExplicitLQHelper( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitLQ( Matrix<F>& L, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitLQ");
#endif
    internal::ExplicitLQHelper( L, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
ExplicitLQ( DistMatrix<F>& L, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ExplicitLQ");
#endif
    internal::ExplicitLQHelper( L, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_EXPLICITLQ_HPP
