/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SVT_NORMAL_HPP
#define ELEM_SVT_NORMAL_HPP

#include ELEM_DIAGONALSCALE_INC
#include ELEM_MAXNORM_INC
#include ELEM_ZERONORM_INC
#include ELEM_SVD_INC
#include ELEM_SOFTTHRESHOLD_INC

namespace elem {

namespace svt {

template<typename F>
inline Int
Normal( Matrix<F>& A, Base<F> tau, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("svt::Normal"))
    typedef Base<F> Real;
    Matrix<F> U( A );
    Matrix<Real> s;
    Matrix<F> V;

    SVD( U, s, V );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

template<typename F>
inline Int
Normal( DistMatrix<F>& A, Base<F> tau, bool relative=false )
{
    DEBUG_ONLY(CallStackEntry cse("svt::Normal"))
    typedef Base<F> Real;
    DistMatrix<F> U( A );
    DistMatrix<Real,VR,STAR> s( A.Grid() );
    DistMatrix<F> V( A.Grid() );

    SVD( U, s, V );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace elem

#endif // ifndef ELEM_SVT_NORMAL_HPP
