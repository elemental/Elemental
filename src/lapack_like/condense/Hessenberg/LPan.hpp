/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESSENBERG_LPAN_HPP
#define EL_HESSENBERG_LPAN_HPP

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
        CSE cse("hessenberg::LPan");
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
        const Range<Int> ind0( 0,   k   ),
                         ind1( k,   k+1 ),
                         ind2( k+1, n   );

        auto a12 = A( ind1, ind2    );
        auto a1  = A( ind1, IR(0,n) );
        auto A2  = A( ind2, IR(0,n) );

        auto alpha12L = A( ind1, IR(k+1,k+2) );
        auto a12R     = A( ind1, IR(k+2,n)   ); 

        // NOTE: Transposes of the horizontal Householder vectors
        auto U0  = U( IR(0,n), ind0 );
        auto u10 = U( ind1,    ind0 );
        auto u21 = U( ind2,    ind1 );
        auto U20 = U( ind2,    ind0 );

        auto V0 = V( IR(0,n), ind0 );
        auto v1 = V( IR(0,n), ind1 );

        auto G00     = G( ind0, ind0 );
        auto g01     = G( ind0, ind1 );
        auto gamma11 = G( ind1, ind1 );

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
        CSE cse("hessenberg::UPan");
        AssertSameGrids( A, t, U_MC_STAR, U_MR_STAR, V_MR_STAR, G_STAR_STAR );
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
        const Range<Int> ind0( 0,   k   ),
                         ind1( k,   k+1 ),
                         ind2( k+1, n   );

        auto a12 = A( ind1, ind2    );
        auto a1  = A( ind1, IR(0,n) );
        auto A2  = A( ind2, IR(0,n) );

        auto alpha12L = A( ind1, IR(k+1,k+2) );
        auto a12R     = A( ind1, IR(k+2,n)   );

        // NOTE: Transposes of the horizontal Householder vectors

        auto U20_MC_STAR = U_MC_STAR( ind2, ind0 );
        auto u21_MC_STAR = U_MC_STAR( ind2, ind1 );

        auto U0_MR_STAR  = U_MR_STAR( IR(0,n), ind0 );
        auto u10_MR_STAR = U_MR_STAR( ind1,    ind0 );
        auto u21_MR_STAR = U_MR_STAR( ind2,    ind1 );

        auto V0_MR_STAR = V_MR_STAR( IR(0,n), ind0 ); 
        auto v1_MR_STAR = V_MR_STAR( IR(0,n), ind1 );

        auto G00_STAR_STAR     = G_STAR_STAR( ind0, ind0 );
        auto g01_STAR_STAR     = G_STAR_STAR( ind0, ind1 );
        auto gamma11_STAR_STAR = G_STAR_STAR( ind1, ind1 );

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
        El::AllReduce( y10_STAR_STAR, U0_MR_STAR.ColComm() );
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
        El::AllReduce( v1_MR_STAR, A2.ColComm() );

        // g01 := U20^H u21
        LocalGemv
        ( ADJOINT, F(1), U20_MC_STAR, u21_MC_STAR, F(0), g01_STAR_STAR );
        El::AllReduce( g01_STAR_STAR, U20_MC_STAR.ColComm() );

        // gamma11 := 1/tau
        gamma11_STAR_STAR.Set(0,0,F(1)/tau);
    }
}

} // namespace hessenberg
} // namespace El

#endif // ifndef EL_HESSENBERG_LPAN_HPP
