/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SVT_CROSS_HPP
#define EL_SVT_CROSS_HPP

namespace El {

namespace svt {

template<typename F>
Int Cross( Matrix<F>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svt::Cross"))
    typedef Base<F> Real;
    Matrix<F> U( A );
    Matrix<Real> s;
    Matrix<F> V;

    svd::Thresholded( U, s, V, tau, relative );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

template<typename F>
Int Cross( DistMatrix<F>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svt::Cross"))
    typedef Base<F> Real;
    DistMatrix<F> U( A );
    DistMatrix<Real,VR,STAR> s( A.Grid() );
    DistMatrix<F> V( A.Grid() );

    svd::Thresholded( U, s, V, tau, relative );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

template<typename F>
Int TallCross( DistMatrix<F,VC,STAR>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("svt::TallCross"))
    typedef Base<F> Real;
    DistMatrix<F,VC,STAR> U( A );
    DistMatrix<Real,STAR,STAR> s( A.Grid() );
    DistMatrix<F,STAR,STAR> V( A.Grid() );

    svd::TallThresholded( U, s, V, tau, relative );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    LocalGemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace El

#endif // ifndef EL_SVT_CROSS_HPP
