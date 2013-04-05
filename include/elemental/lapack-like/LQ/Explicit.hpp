/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_LQ_EXPLICIT_HPP
#define LAPACK_LQ_EXPLICIT_HPP

#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/LQ.hpp"
#include "elemental/matrices/Identity.hpp"

namespace elem {
namespace lq {

template<typename Real>
inline void
ExplicitHelper( Matrix<Real>& A )
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
ExplicitHelper( DistMatrix<Real>& A )
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
ExplicitHelper( Matrix<Complex<Real> >& A )
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
ExplicitHelper( DistMatrix<Complex<Real> >& A )
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
ExplicitHelper( Matrix<Real>& L, Matrix<Real>& A )
{
    LQ( A );
    L = A;
    MakeTriangular( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Real> Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitHelper( DistMatrix<Real>& L, DistMatrix<Real>& A )
{
    LQ( A );
    L = A;
    MakeTriangular( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    const Grid& g = A.Grid();
    DistMatrix<Real> Q( g );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitHelper( Matrix<Complex<Real> >& L, Matrix<Complex<Real> >& A )
{
    Matrix<Complex<Real> > t;
    LQ( A, t );
    L = A;
    MakeTriangular( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    Matrix<Complex<Real> > Q;
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename Real>
inline void
ExplicitHelper
( DistMatrix<Complex<Real> >& L, DistMatrix<Complex<Real> >& A )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MD,STAR> t( g );
    LQ( A, t );
    L = A;
    MakeTriangular( LOWER, L );

    // TODO: Replace this with an in-place expansion of Q
    DistMatrix<Complex<Real> > Q( g );
    Identity( A.Height(), A.Width(), Q );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Q );
    A = Q;
}

template<typename F> 
inline void
Explicit( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("lq::Explicit");
#endif
    ExplicitHelper( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
Explicit( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("lq::Explicit");
#endif
    ExplicitHelper( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
Explicit( Matrix<F>& L, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("lq::Explicit");
#endif
    ExplicitHelper( L, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
Explicit( DistMatrix<F>& L, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("lq::Explicit");
#endif
    ExplicitHelper( L, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace lq
} // namespace elem

#endif // ifndef LAPACK_LQ_EXPLICIT_HPP
