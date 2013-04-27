/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP
#define CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/lapack-like/ApplyRowPivots.hpp"
#include "elemental/lapack-like/QR/BusingerGolub.hpp"
#include "elemental/lapack-like/SVD.hpp"
#include "elemental/convex/SoftThreshold.hpp"

namespace elem {

template<typename F>
inline int
SingularValueSoftThreshold( Matrix<F>& A, BASE(F) tau )
{
#ifndef RELEASE
    CallStackEntry entry("SingularValueSoftThreshold");
#endif
    typedef BASE(F) R;
    Matrix<F> U( A );
    Matrix<R> s;
    Matrix<F> V;

    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

// Preprocess with numSteps iterations of pivoted QR factorization
template<typename Real>
inline int
SingularValueSoftThreshold( Matrix<Real>& A, Real tau, int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("SingularValueSoftThreshold");
    if( numSteps > std::min(A.Height(),A.Width()) )
        throw std::logic_error("number of steps is too large");
#endif
    const int m = A.Height();
    const int n = A.Width();
    Matrix<Real> ACopy( A );
    Matrix<int> p;
    qr::BusingerGolub( ACopy, p, numSteps );
    Matrix<Real> ACopyUpper;
    LockedView( ACopyUpper, ACopy, 0, 0, numSteps, n );

    Matrix<Real> U( ACopyUpper ), s, V;
    MakeTriangular( UPPER, U );
    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    ApplyInverseRowPivots( V, p );
    Matrix<Real> RThresh;
    Gemm( NORMAL, ADJOINT, Real(1), U, V, RThresh );

    ACopy.ResizeTo( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, ACopy );
    Gemm( NORMAL, NORMAL, Real(1), ACopy, RThresh, Real(0), A );

    return ZeroNorm( s );
}

template<typename Real>
inline int
SingularValueSoftThreshold
( Matrix<Complex<Real> >& A, Real tau, int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("SingularValueSoftThreshold");
    if( numSteps > std::min(A.Height(),A.Width()) )
        throw std::logic_error("number of steps is too large");
#endif
    typedef Complex<Real> C;
    const int m = A.Height();
    const int n = A.Width();
    Matrix<C> ACopy( A ), t;
    Matrix<int> p;
    qr::BusingerGolub( ACopy, t, p, numSteps );
    Matrix<C> ACopyUpper;
    LockedView( ACopyUpper, ACopy, 0, 0, numSteps, n );

    Matrix<C> U( ACopyUpper ), V;
    Matrix<Real> s;
    MakeTriangular( UPPER, U );
    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    ApplyInverseRowPivots( V, p );
    Matrix<C> RThresh;
    Gemm( NORMAL, ADJOINT, C(1), U, V, RThresh );

    ACopy.ResizeTo( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, ACopy, t );
    Gemm( NORMAL, NORMAL, C(1), ACopy, RThresh, C(0), A );

    return ZeroNorm( s );
}

template<typename F>
inline int
SingularValueSoftThreshold( DistMatrix<F>& A, BASE(F) tau )
{
#ifndef RELEASE
    CallStackEntry entry("SingularValueSoftThreshold");
#endif
    typedef BASE(F) R;
    DistMatrix<F> U( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    DistMatrix<F> V( A.Grid() );

    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

// Preprocess with numSteps iterations of pivoted QR factorization
template<typename Real>
inline int
SingularValueSoftThreshold( DistMatrix<Real>& A, Real tau, int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("SingularValueSoftThreshold");
    if( numSteps > std::min(A.Height(),A.Width()) )
        throw std::logic_error("number of steps is too large");
#endif
    const int m = A.Height();
    const int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<Real> ACopy( A );
    DistMatrix<int,VR,STAR> p(g);
    qr::BusingerGolub( ACopy, p, numSteps );
    DistMatrix<Real> ACopyUpper(g);
    LockedView( ACopyUpper, ACopy, 0, 0, numSteps, n );

    DistMatrix<Real> U( ACopyUpper ), V(g);
    DistMatrix<Real,VR,STAR> s(g);
    MakeTriangular( UPPER, U );
    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    ApplyInverseRowPivots( V, p );
    DistMatrix<Real> RThresh(g);
    Gemm( NORMAL, ADJOINT, Real(1), U, V, RThresh );

    ACopy.ResizeTo( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, 0, ACopy );
    Gemm( NORMAL, NORMAL, Real(1), ACopy, RThresh, Real(0), A );

    return ZeroNorm( s );
}

template<typename Real>
inline int
SingularValueSoftThreshold
( DistMatrix<Complex<Real> >& A, Real tau, int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("SingularValueSoftThreshold");
    if( numSteps > std::min(A.Height(),A.Width()) )
        throw std::logic_error("number of steps is too large");
#endif
    typedef Complex<Real> C;
    const int m = A.Height();
    const int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<C> ACopy( A );
    DistMatrix<C,MD,STAR> t(g);
    DistMatrix<int,VR,STAR> p(g);
    qr::BusingerGolub( ACopy, t, p, numSteps );
    DistMatrix<C> ACopyUpper(g);
    LockedView( ACopyUpper, ACopy, 0, 0, numSteps, n );

    DistMatrix<C> U( ACopyUpper ), V(g);
    DistMatrix<Real,VR,STAR> s(g);
    MakeTriangular( UPPER, U );
    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    ApplyInverseRowPivots( V, p );
    DistMatrix<C> RThresh(g);
    Gemm( NORMAL, ADJOINT, C(1), U, V, RThresh );

    ACopy.ResizeTo( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, ACopy, t );
    Gemm( NORMAL, NORMAL, C(1), ACopy, RThresh, C(0), A );

    return ZeroNorm( s );
}

} // namespace elem

#endif // ifndef CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP
