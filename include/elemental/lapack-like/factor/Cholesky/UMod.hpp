/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CHOLESKY_UMOD_HPP
#define ELEM_CHOLESKY_UMOD_HPP

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
UUpdate( Matrix<F>& U, Matrix<F>& V )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::mod::UUpdate");
        if( U.Height() != U.Width() )
            LogicError("Cholesky factors must be square");
        if( V.Height() != U.Height() )
            LogicError("V is the wrong height");
    )
    typedef Base<F> Real;
    const Int m = V.Height();
    const Int n = V.Width();

    Matrix<F> z12;

    F* UBuf = U.Buffer();
    const Int ldu = U.LDim();
    for( Int k=0; k<m; ++k )
    {
        F& upsilon11 = UBuf[k+k*ldu];
        auto u12 = ViewRange( U, k, k+1, k+1, m );

        auto v1 = ViewRange( V, k,   0, k+1, n );
        auto V2 = ViewRange( V, k+1, 0, m,   n );

        // Find tau and w such that
        //  /I - tau | 1 | | 1 w^H | \ | upsilon11 | = | beta |
        //  \        | w |           / | v1^H      |   | 0    |
        Conjugate( v1 );
        const F tau = LeftReflector( upsilon11, v1 );

        // Apply the Householder reflector from the left:
        // | u12  | := | u12  | - tau | 1 | | 1 w^H | | u12  |
        // | V2^H |    | V2^H |       | w |           | V2^H |
        //
        //           = | u12  | - tau | 1 | (u12 + w^H V2^H)
        //             | V2^H |       | w |
        // Thus,
        //    u12 -= tau (u12 + w^H V2^H), and
        //    V2^H -= tau w (u12 + w^H V2^H),
        //    V2   -= conj(tau) (u12^H + V2 w) w^H
        //               
        Conjugate( u12, z12 );
        Gemv( NORMAL, F(1), V2, v1, F(1), z12 );
        Ger( -Conj(tau), z12, v1, V2 );
        Conjugate( z12 );
        Axpy( -tau, z12, u12 );
    }
}

template<typename F>
inline void
UDowndate( Matrix<F>& U, Matrix<F>& V )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::mod::UDowndate");
        if( U.Height() != U.Width() )
            LogicError("Cholesky factors must be square");
        if( V.Height() != U.Height() )
            LogicError("V is the wrong height");
    )
    typedef Base<F> Real;
    const Int m = V.Height();
    const Int n = V.Width();

    Matrix<F> z12;

    F* UBuf = U.Buffer();
    const Int ldu = U.LDim();
    for( Int k=0; k<m; ++k )
    {
        F& upsilon11 = UBuf[k+k*ldu];
        auto u12 = ViewRange( U, k, k+1, k+1, m );

        auto v1 = ViewRange( V, k,   0, k+1, n );
        auto V2 = ViewRange( V, k+1, 0, m,   n );

        // Find tau and w such that
        //  /I - 1/tau | 1 | | 1 w^H | Sigma \ | upsilon11 | = | -beta |
        //  \          | w |                 / | v1^H      |   | 0     |
        // where Sigma = diag(+1,-1,...,-1) and beta >= 0
        Conjugate( v1 );
        const F tau = LeftHyperbolicReflector( upsilon11, v1 );

        // Apply the negative of the Householder reflector from the left:
        // | u12  | := -| u12  | + 1/tau | 1 | | 1 w^H | Sigma | u12  |
        // | V2^H |     | V2^H |         | w |                 | V2^H |
        //
        //           = -| u12  | + 1/tau | 1 | (u12 - w^H V2^H)
        //              | V2^H |         | w |
        // Thus,
        //    u12 := -u12 + 1/tau (u12 - w^H V2^H), and
        //    V2^H := -V2^H + 1/tau w (u12 - w^H V2^H),
        //    V2   := -V2 + 1/tau (u12^H - V2 w) w^H
        //               
        upsilon11 = -upsilon11;
        Conjugate( u12, z12 );
        Gemv( NORMAL, F(-1), V2, v1, F(1), z12 );
        Scale( F(-1), u12 );
        Scale( F(-1), V2 );
        Ger( F(1)/tau, z12, v1, V2 );
        Conjugate( z12 );
        Axpy( F(1)/tau, z12, u12 );
    }
}

} // namespace mod

template<typename F>
inline void
UMod( Matrix<F>& U, Base<F> alpha, Matrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("cholesky::UMod"))
    typedef Base<F> Real;
    if( alpha == Real(0) )
        return;
    else if( alpha > Real(0) )
    {
        Scale( Sqrt(alpha), V );
        mod::UUpdate( U, V );
    }
    else
    {
        Scale( Sqrt(-alpha), V );
        mod::UDowndate( U, V );
    }
}

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_CHOLESKY_UMOD_HPP
