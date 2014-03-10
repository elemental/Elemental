/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CHOLESKY_LMOD_HPP
#define ELEM_CHOLESKY_LMOD_HPP

#include ELEM_AXPY_INC
#include ELEM_GEMV_INC
#include ELEM_GER_INC
#include ELEM_REFLECTOR_INC
#include ELEM_HYPERBOLICREFLECTOR_INC

// TODO: Blocked and distributed algorithms

namespace elem {
namespace cholesky {

namespace mod {

template<typename F>
inline void
LUpdate( Matrix<F>& L, Matrix<F>& V )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::mod::LUpdate");
        if( L.Height() != L.Width() )
            LogicError("Cholesky factors must be square");
        if( V.Height() != L.Height() )
            LogicError("V is the wrong height");
    )
    typedef Base<F> Real;
    const Int m = V.Height();
    const Int n = V.Width();

    Matrix<F> z21;

    F* LBuf = L.Buffer();
    const Int ldl = L.LDim();
    for( Int k=0; k<m; ++k )
    {
        F& lambda11 = LBuf[k+k*ldl];
        auto l21 = ViewRange( L, k+1, k, m, k+1 );

        auto v1 = View( V, k, 0, 1, n );
        auto V2 = ViewRange( V, k+1, 0, m, n );

        // Find tau and u such that
        //  | lambda11 u | /I - tau | 1   | | 1 conj(u) |\ = | beta 0 |
        //                 \        | u^T |              /
        const F tau = RightReflector( lambda11, v1 );

        // Apply the Householder reflector from the right:
        // | l21 V2 | := | l21 V2 | - tau | l21 V2 | | 1   | | 1 conj(u) |
        //                                           | u^T |
        //             = | l21 V2 | - tau (l21 + V2 u^T) | 1 conj(u) | 
        z21 = l21;
        Gemv( NORMAL, F(1), V2, v1, F(1), z21 );
        Axpy( -tau, z21, l21 );
        Ger( -tau, z21, v1, V2 );
    }
}

template<typename F>
inline void
LDowndate( Matrix<F>& L, Matrix<F>& V )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::mod::LDowndate");
        if( L.Height() != L.Width() )
            LogicError("Cholesky factors must be square");
        if( V.Height() != L.Height() )
            LogicError("V is the wrong height");
    )
    typedef Base<F> Real;
    const Int m = V.Height();
    const Int n = V.Width();

    Matrix<F> z21;

    F* LBuf = L.Buffer();
    const Int ldl = L.LDim();
    for( Int k=0; k<m; ++k )
    {
        F& lambda11 = LBuf[k+k*ldl];
        auto l21 = ViewRange( L, k+1, k, m, k+1 );

        auto v1 = View( V, k, 0, 1, n );
        auto V2 = ViewRange( V, k+1, 0, m, n );

        // Find tau and u such that
        //  | lambda11 u | /I - 1/tau Sigma | 1   | | 1 conj(u) |\ = | -beta 0 |
        //                 \                | u^T |              /
        // where Sigma = diag(+1,-1,...,-1) and beta >= 0
        const F tau = RightHyperbolicReflector( lambda11, v1 );

        // Apply the negative of the hyperbolic Householder reflector from the 
        // right:
        // |l21 V2| := -|l21 V2| + 1/tau |l21 V2| Sigma |1  | |1 conj(u)|
        //                                              |u^T|
        //           = -|l21 V2| + 1/tau |l21 -V2| |1  | |1 conj(u)|
        //                                         |u^T|
        //           = -|l21 V2| + 1/tau |l21 - V2 u^T| |1 conj(u)|
        lambda11 = -lambda11;
        z21 = l21;
        Gemv( NORMAL, F(-1), V2, v1, F(1), z21 );
        Scale( F(-1), l21 );
        Scale( F(-1), V2  );
        Axpy( F(1)/tau, z21, l21 );
        Ger( F(1)/tau, z21, v1, V2 );
    }
}

} // namespace mod

template<typename F>
inline void
LMod( Matrix<F>& L, Base<F> alpha, Matrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::LMod"))
    typedef Base<F> Real;
    if( alpha == Real(0) )
        return;
    else if( alpha > Real(0) )
    {
        Scale( Sqrt(alpha), V );  
        mod::LUpdate( L, V );
    }
    else
    {
        Scale( Sqrt(-alpha), V );
        mod::LDowndate( L, V );
    }
}

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_CHOLESKY_LMOD_HPP
