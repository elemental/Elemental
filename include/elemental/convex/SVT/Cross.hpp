/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONVEX_SVT_CROSS_HPP
#define ELEM_CONVEX_SVT_CROSS_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/lapack-like/Norm/Zero.hpp"
#include "elemental/lapack-like/SVD.hpp"
#include "elemental/convex/SoftThreshold.hpp"

namespace elem {

namespace svt {

template<typename F>
inline Int
Cross( Matrix<F>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry entry("svt::Cross");
#endif
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
inline Int
Cross( DistMatrix<F>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry entry("svt::Cross");
#endif
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
inline Int
TallCross( DistMatrix<F,VC,STAR>& A, BASE(F) tau, bool relative=false )
{
#ifndef RELEASE
    CallStackEntry entry("svt::TallCross");
#endif
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
} // namespace elem

#endif // ifndef ELEM_CONVEX_SVT_CROSS_HPP
