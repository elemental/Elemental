/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HESSENBERG_UPAN_HPP
#define ELEM_HESSENBERG_UPAN_HPP

#include ELEM_GEMV_INC
#include ELEM_TRSV_INC

#include ELEM_REFLECTOR_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace hessenberg {

// NOTE: This is an extension into complex arithmetic of the sequential algorithm proposed in:
//       G. Quintana-Orti and R. van de Geijn,  
//       "Improving the performance of reduction to Hessenberg form" 
//       After switching to complex arithmetic, it was more natural to switch 
//       to lower-triangular matrices in the UT transform.

// NOTE: It would be possible to avoid the space for U if we were more careful
//       about applying the portion interleaved with the Hessenberg matrix.
template<typename F>
inline void UPan
( Matrix<F>& A, Matrix<F>& t, Matrix<F>& U, Matrix<F>& V, Matrix<F>& G )
{
    const Int nU = U.Width();
    const Int n = A.Height();
    DEBUG_ONLY(
        CallStackEntry cse("hessenberg::UPan");
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
        Gemv( NORMAL, F(1), A2, u21, F(0), v1 );

        // g10 := u21^H U20 = (U20^H u21)^H
        Gemv( ADJOINT, F(1), U20, u21, F(0), g10 );
        Conjugate( g10 );
        
        // gamma11 := 1/tau
        gamma11.Set(0,0,F(1)/tau);
    }
}

template<typename F>
inline void UPan
( DistMatrix<F>& A, 
  DistMatrix<F,STAR,STAR>& t, 
  DistMatrix<F,MC,  STAR>& U_MC_STAR, 
  DistMatrix<F,MR,  STAR>& U_MR_STAR,
  DistMatrix<F,MC,  STAR>& V_MC_STAR, 
  DistMatrix<F,STAR,STAR>& G_STAR_STAR )
{
    const Int nU = U_MC_STAR.Width();
    const Int n = A.Height();
    DEBUG_ONLY(
        CallStackEntry cse("hessenberg::UPan");
        if( A.Grid() != t.Grid() || t.Grid() != U_MC_STAR.Grid() || 
            U_MC_STAR.Grid() != U_MR_STAR.Grid() || 
            U_MR_STAR.Grid() != V_MC_STAR.Grid() ||
            V_MC_STAR.Grid() != G_STAR_STAR.Grid() )
            LogicError("Grids must match");
        if( A.ColAlign() != U_MC_STAR.ColAlign() )
            LogicError("A and U[MC,* ] must be aligned");
        if( A.RowAlign() != U_MR_STAR.ColAlign() )
            LogicError("A and U[MR,* ] must be aligned");
        if( A.ColAlign() != V_MC_STAR.ColAlign() )
            LogicError("A and V[MC,* ] must be aligned");
        if( nU >= n )            
            LogicError("V is too wide for the panel factorization");
        if( U_MC_STAR.Height() != A.Height() )
            LogicError("U[MC,* ] must be the same height as A");
        if( U_MR_STAR.Height() != A.Height() )
            LogicError("U[MR,* ] must be the same height as A");
        if( U_MR_STAR.Width() != nU )
            LogicError("U[MR,* ] must be the same width as U[MC,* ]");
        if( V_MC_STAR.Height() != A.Height() )
            LogicError("V[MC,* ] must be the same height as A");
        if( V_MC_STAR.Width() != nU )
            LogicError("V[MC,* ] must be the same width as U");
    )
    const Grid& g = A.Grid();

    Zeros( U_MC_STAR,   n,  nU );
    Zeros( U_MR_STAR,   n,  nU );
    Zeros( V_MC_STAR,   n,  nU );
    Zeros( G_STAR_STAR, nU, nU );

    DistMatrix<F,MC,  STAR> a1_MC_STAR(g);
    DistMatrix<F,STAR,STAR> y10_STAR_STAR(g);

    for( Int k=0; k<nU; ++k )
    {
        auto a21      = ViewRange( A, k+1, k,   n,   k+1 );
        auto alpha21T = ViewRange( A, k+1, k,   k+2, k+1 );
        auto a21B     = ViewRange( A, k+2, k,   n,   k+1 );
        auto a1       = ViewRange( A, 0,   k,   n,   k+1 );
        auto A2       = ViewRange( A, 0,   k+1, n,   n   );

        auto U0_MC_STAR  = ViewRange( U_MC_STAR, 0,   0, n,   k   );
        auto u10_MC_STAR = ViewRange( U_MC_STAR, k,   0, k+1, k   );
        auto u21_MC_STAR = ViewRange( U_MC_STAR, k+1, k, n,   k+1 );
        auto u21_MR_STAR = ViewRange( U_MR_STAR, k+1, k, n,   k+1 );
        auto U20_MR_STAR = ViewRange( U_MR_STAR, k+1, 0, n,   k   );

        auto V0_MC_STAR = ViewRange( V_MC_STAR, 0, 0, n, k   );
        auto v1_MC_STAR = ViewRange( V_MC_STAR, 0, k, n, k+1 );

        auto G00_STAR_STAR     = ViewRange( G_STAR_STAR, 0, 0, k,   k   );
        auto g10_STAR_STAR     = ViewRange( G_STAR_STAR, k, 0, k+1, k   );
        auto gamma11_STAR_STAR = ViewRange( G_STAR_STAR, k, k, k+1, k+1 );

        // a1 := (I - U0 inv(G00) U0^H) (a1 - V0 inv(G00)^H u10^H)
        // -------------------------------------------------------
        // a1 := a1 - V0 inv(G00)^H u10^H
        a1_MC_STAR.AlignWith( a1 );
        a1_MC_STAR = a1;
        Conjugate( u10_MC_STAR, y10_STAR_STAR );
        Trsv
        ( LOWER, ADJOINT, NON_UNIT, 
          G00_STAR_STAR.LockedMatrix(), y10_STAR_STAR.Matrix() );
        LocalGemv( NORMAL, F(-1), V0_MC_STAR, y10_STAR_STAR, F(1), a1_MC_STAR );
        // a1 := a1 - U0 (inv(G00) (U0^H a1))
        LocalGemv
        ( ADJOINT, F(1), U0_MC_STAR, a1_MC_STAR, F(0), y10_STAR_STAR ); 
        y10_STAR_STAR.SumOver( U0_MC_STAR.ColComm() );
        Trsv
        ( LOWER, NORMAL, NON_UNIT, 
          G00_STAR_STAR.LockedMatrix(), y10_STAR_STAR.Matrix() );
        LocalGemv( NORMAL, F(-1), U0_MC_STAR, y10_STAR_STAR, F(1), a1_MC_STAR );
        a1 = a1_MC_STAR; 

        // Find tau and v such that
        //  / I - tau | 1 | | 1, v^H | \ | alpha21T | = | beta |
        //  \         | v |            / |     a21B |   |    0 |
        const F tau = LeftReflector( alpha21T, a21B );
        t.Set(k,0,tau);

        // Store u21 := | 1 |
        //              | v |
        u21_MC_STAR = a21;
        u21_MR_STAR = a21;
        u21_MC_STAR.Set(0,0,F(1));
        u21_MR_STAR.Set(0,0,F(1));

        // v1 := A2 u21
        LocalGemv( NORMAL, F(1), A2, u21_MR_STAR, F(0), v1_MC_STAR );
        v1_MC_STAR.SumOver( A2.RowComm() );

        // g10 := u21^H U20 = (U20^H u21)^H
        LocalGemv
        ( ADJOINT, F(1), U20_MR_STAR, u21_MR_STAR, F(0), g10_STAR_STAR );
        g10_STAR_STAR.SumOver( U20_MR_STAR.ColComm() );
        Conjugate( g10_STAR_STAR );
        
        // gamma11 := 1/tau
        gamma11_STAR_STAR.Set(0,0,F(1)/tau);
    }
}

} // namespace hessenberg
} // namespace elem

#endif // ifndef ELEM_HESSENBERG_UPAN_HPP
