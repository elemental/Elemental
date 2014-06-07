/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HESSENBERG_LPAN_HPP
#define EL_HESSENBERG_LPAN_HPP

#include EL_ZEROS_INC

namespace El {
namespace hessenberg {

// NOTE: This is an extension into complex arithmetic 
//       (and conjugate-transposition) of the sequential algorithm proposed in:
//       G. Quintana-Orti and R. van de Geijn,  
//       "Improving the performance of reduction to Hessenberg form" 

// NOTE: It would be possible to avoid the space for U if we were more careful
//       about applying the portion interleaved with the Hessenberg matrix.
template<typename F>
inline void LPan
( Matrix<F>& A, Matrix<F>& t, Matrix<F>& U, Matrix<F>& V, Matrix<F>& G )
{
    const Int nU = U.Width();
    const Int n = A.Height();
    DEBUG_ONLY(
        CallStackEntry cse("hessenberg::LPan");
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
        auto a12      = ViewRange( A, k,   k+1, k+1, n   );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+2 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, n   );
        auto a1       = ViewRange( A, k,   0,   k+1, n   );
        auto A2       = ViewRange( A, k+1, 0,   n,   n   ); 

        // NOTE: Transposes of the horizontal Householder vectors
        auto U0  = ViewRange( U, 0,   0, n,   k   );
        auto u10 = ViewRange( U, k,   0, k+1, k );
        auto u21 = ViewRange( U, k+1, k, n,   k+1 );
        auto U20 = ViewRange( U, k+1, 0, n,   k   );

        auto V0 = ViewRange( V, 0, 0, n, k   );
        auto v1 = ViewRange( V, 0, k, n, k+1 );

        auto G00     = ViewRange( G, 0, 0, k, k   );
        auto g01     = ViewRange( G, 0, k, k, k+1 );
        auto gamma11 = ViewRange( G, k, k, k+1, k+1 );

        // a1 := (a1 - u10 inv(G00)^H V0^H) (I - U0 inv(G00) U0^H)
        // -------------------------------------------------------
        // a1 := conj(a1 - (u10 inv(G00)^H) V0^H)
        //     = conj(a1 - (V0 (inv(G00) u10^H))^H)
        //     = conj(a1) - (V0 (inv(G00) u10^H))^T
        Conjugate( u10, y10 ); 
        Trsv( UPPER, NORMAL, NON_UNIT, G00, y10 );
        Conjugate( a1 );
        Gemv( NORMAL, F(-1), V0, y10, F(1), a1 );
        // a1 := conj(a1) - conj(a1) U0 inv(G00) U0^H
        //     = conj(a1 - (U0 (inv(G00)^H (U0^H a1^T)))^T)
        Gemv( ADJOINT, F(1), U0, a1, F(0), y10 );
        Trsv( UPPER, ADJOINT, NON_UNIT, G00, y10 );
        Gemv( NORMAL, F(-1), U0, y10, F(1), a1 );
        Conjugate( a1 );

        // Find tau and v such that
        // | alpha12L a12R | / I - tau | 1   | | 1 conj(v) | \ = | beta 0 |
        //                   \         | v^T |               /
        const F tau = RightReflector( alpha12L, a12R );
        t.Set(k,0,tau);

        // Store u21 := | 1   |
        //              | v^T |
        Transpose( a12, u21 );
        u21.Set(0,0,F(1));

        // v1 := A2^H u21
        Gemv( ADJOINT, F(1), A2, u21, F(0), v1 );

        // g01 := U20^H u21
        Gemv( ADJOINT, F(1), U20, u21, F(0), g01 );
        
        // gamma11 := 1/tau
        gamma11.Set(0,0,F(1)/tau);
    }
}

template<typename F>
inline void LPan
( DistMatrix<F>& A,
  DistMatrix<F,STAR,STAR>& t,
  DistMatrix<F,MC,  STAR>& U_MC_STAR,
  DistMatrix<F,MR,  STAR>& U_MR_STAR,
  DistMatrix<F,MR,  STAR>& V_MR_STAR,
  DistMatrix<F,STAR,STAR>& G_STAR_STAR )
{
    const Int nU = U_MC_STAR.Width();
    const Int n = A.Height();
    DEBUG_ONLY(
        CallStackEntry cse("hessenberg::UPan");
        if( A.Grid() != t.Grid() || t.Grid() != U_MC_STAR.Grid() ||
            U_MC_STAR.Grid() != U_MR_STAR.Grid() ||
            U_MR_STAR.Grid() != V_MR_STAR.Grid() ||
            V_MR_STAR.Grid() != G_STAR_STAR.Grid() )
            LogicError("Grids must match");
        if( A.ColAlign() != U_MC_STAR.ColAlign() )
            LogicError("A and U[MC,* ] must be aligned");
        if( A.RowAlign() != U_MR_STAR.ColAlign() )
            LogicError("A and U[MR,* ] must be aligned");
        if( A.RowAlign() != V_MR_STAR.ColAlign() )
            LogicError("A and V[MR,* ] must be aligned");
        if( nU >= n )
            LogicError("V is too wide for the panel factorization");
        if( U_MC_STAR.Height() != A.Height() )
            LogicError("U[MC,* ] must be the same height as A");
        if( U_MR_STAR.Height() != A.Height() )
            LogicError("U[MR,* ] must be the same height as A");
        if( U_MR_STAR.Width() != nU )
            LogicError("U[MR,* ] must be the same width as U[MC,* ]");
        if( V_MR_STAR.Height() != A.Height() )
            LogicError("V[MR,* ] must be the same height as A");
        if( V_MR_STAR.Width() != nU )
            LogicError("V[MR,* ] must be the same width as U");
    )
    const Grid& g = A.Grid();

    Zeros( U_MC_STAR,   n,  nU );
    Zeros( U_MR_STAR,   n,  nU );
    Zeros( V_MR_STAR,   n,  nU );
    Zeros( G_STAR_STAR, nU, nU );

    DistMatrix<F,STAR,MR  > a1Conj_STAR_MR(g);
    DistMatrix<F,STAR,STAR> y10_STAR_STAR(g);

    for( Int k=0; k<nU; ++k )
    {
        auto a12      = ViewRange( A, k,   k+1, k+1, n   );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+2 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, n   );
        auto a1       = ViewRange( A, k,   0,   k+1, n   );
        auto A2       = ViewRange( A, k+1, 0,   n,   n   );

        // NOTE: Transposes of the horizontal Householder vectors
        auto U0_MR_STAR  = ViewRange( U_MR_STAR, 0,   0, n,   k   );
        auto u10_MR_STAR = ViewRange( U_MR_STAR, k,   0, k+1, k   );
        auto u21_MC_STAR = ViewRange( U_MC_STAR, k+1, k, n,   k+1 );
        auto u21_MR_STAR = ViewRange( U_MR_STAR, k+1, k, n,   k+1 );
        auto U20_MC_STAR = ViewRange( U_MC_STAR, k+1, 0, n,   k   );

        auto V0_MR_STAR = ViewRange( V_MR_STAR, 0, 0, n, k   );
        auto v1_MR_STAR = ViewRange( V_MR_STAR, 0, k, n, k+1 );

        auto G00_STAR_STAR     = ViewRange( G_STAR_STAR, 0, 0, k, k   );
        auto g01_STAR_STAR     = ViewRange( G_STAR_STAR, 0, k, k, k+1 );
        auto gamma11_STAR_STAR = ViewRange( G_STAR_STAR, k, k, k+1, k+1 );

        // a1 := (a1 - u10 inv(G00)^H V0^H) (I - U0 inv(G00) U0^H)
        // -------------------------------------------------------
        // a1 := conj(a1 - (u10 inv(G00)^H) V0^H)
        //     = conj(a1 - (V0 (inv(G00) u10^H))^H)
        //     = conj(a1) - (V0 (inv(G00) u10^H))^T
        a1Conj_STAR_MR.AlignWith( a1 );
        Conjugate( a1, a1Conj_STAR_MR );
        Conjugate( u10_MR_STAR, y10_STAR_STAR );
        Trsv
        ( UPPER, NORMAL, NON_UNIT, 
          G00_STAR_STAR.LockedMatrix(), y10_STAR_STAR.Matrix() );
        LocalGemv
        ( NORMAL, F(-1), V0_MR_STAR, y10_STAR_STAR, F(1), a1Conj_STAR_MR );
        // a1 := conj(a1) - conj(a1) U0 inv(G00) U0^H
        //     = conj(a1 - (U0 (inv(G00)^H (U0^H a1^T)))^T)
        LocalGemv
        ( ADJOINT, F(1), U0_MR_STAR, a1Conj_STAR_MR, F(0), y10_STAR_STAR );
        y10_STAR_STAR.SumOver( U0_MR_STAR.ColComm() );
        Trsv
        ( UPPER, ADJOINT, NON_UNIT, 
          G00_STAR_STAR.LockedMatrix(), y10_STAR_STAR.Matrix() );
        LocalGemv
        ( NORMAL, F(-1), U0_MR_STAR, y10_STAR_STAR, F(1), a1Conj_STAR_MR );
        Conjugate( a1Conj_STAR_MR, a1 ); 

        // Find tau and v such that
        // | alpha12L a12R | / I - tau | 1   | | 1 conj(v) | \ = | beta 0 |
        //                   \         | v^T |               /
        const F tau = RightReflector( alpha12L, a12R );
        t.Set(k,0,tau);

        // Store u21 := | 1   |
        //              | v^T |
        Transpose( a12, u21_MC_STAR );
        Transpose( a12, u21_MR_STAR );
        u21_MC_STAR.Set(0,0,F(1));
        u21_MR_STAR.Set(0,0,F(1));

        // v1 := A2^H u21
        LocalGemv( ADJOINT, F(1), A2, u21_MC_STAR, F(0), v1_MR_STAR );
        v1_MR_STAR.SumOver( A2.ColComm() );

        // g01 := U20^H u21
        LocalGemv
        ( ADJOINT, F(1), U20_MC_STAR, u21_MC_STAR, F(0), g01_STAR_STAR );
        g01_STAR_STAR.SumOver( U20_MC_STAR.ColComm() );

        // gamma11 := 1/tau
        gamma11_STAR_STAR.Set(0,0,F(1)/tau);
    }
}

} // namespace hessenberg
} // namespace El

#endif // ifndef EL_HESSENBERG_LPAN_HPP
