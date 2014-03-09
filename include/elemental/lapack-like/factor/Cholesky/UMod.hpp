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

// TODO: Downdating, blocked algorithms, and distributed algorithms

namespace elem {
namespace cholesky {

template<typename F>
inline void
UMod( Matrix<F>& U, Base<F> alpha, Matrix<F>& V )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UMod");
        if( U.Height() != U.Width() )
            LogicError("Cholesky factors must be square");
        if( V.Height() != U.Height() )
            LogicError("V is the wrong height");
        if( alpha < Base<F>(0) )
            LogicError("Only updates are currently supported");
    )
    typedef Base<F> Real;
    const Int m = V.Height();
    const Int n = V.Width();

    Matrix<F> z12;

    F* UBuf = U.Buffer();
    const Int ldu = U.LDim();
    Scale( Sqrt(alpha), V );
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

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_CHOLESKY_UMOD_HPP
