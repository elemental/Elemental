/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONVEX_SVT_PIVOTEDQR_HPP
#define ELEM_CONVEX_SVT_PIVOTEDQR_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/lapack-like/Norm/Zero.hpp"
#include "elemental/lapack-like/QR.hpp"
#include "elemental/lapack-like/SVD.hpp"
#include "elemental/convex/SoftThreshold.hpp"

namespace elem {
namespace svt {

// Preprocess with numSteps iterations of pivoted QR factorization

template<typename F>
inline Int
PivotedQR( Matrix<F>& A, BASE(F) tau, Int numSteps, bool relative=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("svt::PivotedQR");
        if( numSteps > std::min(A.Height(),A.Width()) )
            LogicError("number of steps is too large");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<F> ACopy( A ), t;
    Matrix<Int> p;
    qr::BusingerGolub( ACopy, t, p, numSteps );
    Matrix<F> ACopyUpper;
    LockedView( ACopyUpper, ACopy, 0, 0, numSteps, n );

    Matrix<F> U( ACopyUpper ), V;
    Matrix<Real> s;
    MakeTriangular( UPPER, U );
    svd::Thresholded( U, s, V, tau, relative );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    ApplyInverseRowPivots( V, p );
    Matrix<F> RThresh;
    Gemm( NORMAL, ADJOINT, F(1), U, V, RThresh );

    ACopy.Resize( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, ACopy, t );
    Gemm( NORMAL, NORMAL, F(1), ACopy, RThresh, F(0), A );

    return ZeroNorm( s );
}


template<typename F>
inline Int
PivotedQR( DistMatrix<F>& A, BASE(F) tau, Int numSteps, bool relative=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("svt::PivotedQR");
        if( numSteps > std::min(A.Height(),A.Width()) )
            LogicError("number of steps is too large");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<F> ACopy( A );
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Int,VR,STAR> p(g);
    qr::BusingerGolub( ACopy, t, p, numSteps );
    DistMatrix<F> ACopyUpper(g);
    LockedView( ACopyUpper, ACopy, 0, 0, numSteps, n );

    DistMatrix<F> U( ACopyUpper ), V(g);
    DistMatrix<Real,VR,STAR> s(g);
    MakeTriangular( UPPER, U );
    svd::Thresholded( U, s, V, tau, relative );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    ApplyInverseRowPivots( V, p );
    DistMatrix<F> RThresh(g);
    Gemm( NORMAL, ADJOINT, F(1), U, V, RThresh );

    ACopy.Resize( m, numSteps );
    ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, ACopy, t );
    Gemm( NORMAL, NORMAL, F(1), ACopy, RThresh, F(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace elem

#endif // ifndef ELEM_CONVEX_SVT_PIVOTEDQR_HPP
