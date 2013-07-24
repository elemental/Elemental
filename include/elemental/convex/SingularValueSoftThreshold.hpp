/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP
#define ELEM_CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP

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
template<typename F>
inline int
SingularValueSoftThreshold( Matrix<F>& A, BASE(F) tau, int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("SingularValueSoftThreshold");
    if( numSteps > std::min(A.Height(),A.Width()) )
        throw std::logic_error("number of steps is too large");
#endif
    typedef BASE(F) Real;
    const int m = A.Height();
    const int n = A.Width();
    Matrix<F> ACopy( A ), t;
    Matrix<int> p;
    qr::BusingerGolub( ACopy, t, p, numSteps );
    Matrix<F> ACopyUpper;
    LockedView( ACopyUpper, ACopy, 0, 0, numSteps, n );

    Matrix<F> U( ACopyUpper ), V;
    Matrix<Real> s;
    MakeTriangular( UPPER, U );
    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    ApplyInverseRowPivots( V, p );
    Matrix<F> RThresh;
    Gemm( NORMAL, ADJOINT, F(1), U, V, RThresh );

    ACopy.ResizeTo( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, ACopy, t );
    Gemm( NORMAL, NORMAL, F(1), ACopy, RThresh, F(0), A );

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
template<typename F>
inline int
SingularValueSoftThreshold( DistMatrix<F>& A, BASE(F) tau, int numSteps )
{
#ifndef RELEASE
    CallStackEntry entry("SingularValueSoftThreshold");
    if( numSteps > std::min(A.Height(),A.Width()) )
        throw std::logic_error("number of steps is too large");
#endif
    typedef BASE(F) Real;
    const int m = A.Height();
    const int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<F> ACopy( A );
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<int,VR,STAR> p(g);
    qr::BusingerGolub( ACopy, t, p, numSteps );
    DistMatrix<F> ACopyUpper(g);
    LockedView( ACopyUpper, ACopy, 0, 0, numSteps, n );

    DistMatrix<F> U( ACopyUpper ), V(g);
    DistMatrix<Real,VR,STAR> s(g);
    MakeTriangular( UPPER, U );
    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    ApplyInverseRowPivots( V, p );
    DistMatrix<F> RThresh(g);
    Gemm( NORMAL, ADJOINT, F(1), U, V, RThresh );

    ACopy.ResizeTo( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, ACopy, t );
    Gemm( NORMAL, NORMAL, F(1), ACopy, RThresh, F(0), A );

    return ZeroNorm( s );
}

} // namespace elem

#endif // ifndef ELEM_CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP
