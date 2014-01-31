/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HESSENBERG_PANELU_HPP
#define ELEM_HESSENBERG_PANELU_HPP

#include ELEM_GEMV_INC
#include ELEM_TRSV_INC

#include ELEM_REFLECTOR_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace hessenberg {

// NOTE: This is an extension into complex arithmetic of the sequential algorithm proposed in:
//       G. Quintana-Orti and R. van de Geijn,  
//       "Improving the performance of reduction to Hessenberg form" 
//       After switching to complex arithmetic, it was more natural to switch to lower-triangular
//       matrices in the UT transform.

// NOTE: It would be possible to avoid the space for U if we were more careful
//       about applying the portion interleaved with the Hessenberg matrix.
template<typename F>
inline void PanelU
( Matrix<F>& A, Matrix<F>& t, Matrix<F>& U, Matrix<F>& V, Matrix<F>& G )
{
    const Int nU = U.Width();
    const Int n = A.Height();
    DEBUG_ONLY(
        CallStackEntry cse("hessenberg::PanelU");
        if( nU >= n )            
            LogicError("V is too wide for the panel factorization");
        if( U.Height() != A.Height() )
            LogicError("U must be the same height as A");
        if( V.Height() != A.Height() )
            LogicError("V must be the same height as A");
        if( V.Width() != nU )
            LogicError("V must be the same width as U");
    )
    const Int tHeight = Max(nU,0);
    t.Resize( tHeight, 1 );

    Zeros( U, n, nU );
    Zeros( V, n, nU );
    Zeros( G, nU, nU );

    Matrix<F> y10;

    for( Int k=0; k<nU; ++k )
    {
        auto a21      = ViewRange( A, k+1, k,   n,   k+1 );
        auto alpha21T = ViewRange( A, k+1, k,   k+2, k+1 );
        auto a21B     = ViewRange( A, k+2, k,   n,   k+1 );
        auto A22      = ViewRange( A, k+1, k+1, n,   n   );
        auto a1       = ViewRange( A, 0,   k,   n,   k+1 );
        auto A2       = ViewRange( A, 0,   k+1, n,   n   );

        auto U0  = ViewRange( U, 0,   0, n,   k   );
        auto u10 = ViewRange( U, k,   0, k+1, k );
        auto u21 = ViewRange( U, k+1, k, n,   k+1 );
        auto U20 = ViewRange( U, k+1, 0, n,   k   );

        auto V0 = ViewRange( V, 0, 0, n, k   );
        auto v1 = ViewRange( V, 0, k, n, k+1 );

        auto G00     = ViewRange( G, 0, 0, k,   k   );
        auto g10     = ViewRange( G, k, 0, k+1, k   );
        auto gamma11 = ViewRange( G, k, k, k+1, k+1 );

        // a1 := (I - U0 inv(G00) U0^H) (a1 - V0 inv(G00)^H u10^H)
        // -------------------------------------------------------
        // a1 := a1 - V0 inv(G00)^H u10^H
        Conjugate( u10, y10 );
        Trsv( LOWER, ADJOINT, NON_UNIT, G00, y10 );
        Gemv( NORMAL, F(-1), V0, y10, F(1), a1 );
        // a1 := a1 - U0 (inv(G00) (U0^H a1))
        Gemv( ADJOINT, F(1), U0, a1, F(0), y10 ); 
        Trsv( LOWER, NORMAL, NON_UNIT, G00, y10 );
        Gemv( NORMAL, F(-1), U0, y10, F(1), a1 );

        // Find tau and v such that
        //  / I - tau | 1 | | 1, v^H | \ | alpha21T | = | beta |
        //  \         | v |            / |     a21B |   |    0 |
        const F tau = LeftReflector( alpha21T, a21B );
        t.Set(k,0,tau);

        // Store u21 := | 1 |
        //              | v |
        u21 = a21;
        u21.Set(0,0,F(1));

        // v1 := A2 u21
        Zeros( v1, n, 1 );
        Gemv( NORMAL, F(1), A2, u21, F(0), v1 );

        // g10 := u21^H U20 = (U20^H u21)^H
        Zeros( g10, 1, k );
        Gemv( ADJOINT, F(1), U20, u21, F(0), g10 );
        Conjugate( g10 );
        
        // gamma11 := 1/tau
        gamma11.Set(0,0,F(1)/tau);
    }
}

} // namespace hessenberg
} // namespace elem

#endif // ifndef ELEM_HESSENBERG_PANELU_HPP
