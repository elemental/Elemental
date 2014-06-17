/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HESSENBERG_LUNB_HPP
#define EL_HESSENBERG_LUNB_HPP

namespace El {
namespace hessenberg {

template<typename F>
inline void LUnb( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::LUnb"))
    const Int n = A.Height();
    const Int tHeight = Max(n-1,0);
    t.Resize( tHeight, 1 );

    Matrix<F> z1, z21;

    for( Int k=0; k<n-1; ++k )
    {
        auto a12      = ViewRange( A, k,   k+1, k+1, n   );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+2 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, n   );
        auto A22      = ViewRange( A, k+1, k+1, n,   n   );
        auto A2       = ViewRange( A, k+1, 0,   n,   n   );

        // Find tau and v such that
        //  |alpha12L a12R| /I - tauP | 1   | | 1 conj(v) |\ = |beta 0|
        //                  \         | v^T |              /
        const F tau = RightReflector( alpha12L, a12R );
        t.Set(k,0,tau);

        // Temporarily set a12 = | 1 v |
        const F beta = alpha12L.Get(0,0);
        alpha12L.Set(0,0,F(1));

        // A2 := Hous(a12^T,tau)^H A2
        //     = (I - conj(tau) a12^T conj(a12)) A2
        //     = A2 - conj(tau) a12^T (A2^H a12^T)^H
        // -----------------------------------------
        // z1 := A2^H a12^T
        Zeros( z1, n, 1 ); 
        Gemv( ADJOINT, F(1), A2, a12, F(0), z1 );
        // A2 := A2 - conj(tau) a12^T z1^H
        Ger( -Conj(tau), a12, z1, A2 );

        // A22 := A22 Hous(a12^T,tau)
        //      = A22 (I - tau a12^T conj(a12))
        //      = A22 - tau (A22 a12^T) conj(a12)
        // --------------------------------------
        // z21 := A22 a12^T
        Zeros( z21, A22.Height(), 1 );
        Gemv( NORMAL, F(1), A22, a12, F(0), z21 );
        // A22 := A22 - tau z21 conj(a12)
        Ger( -tau, z21, a12, A22 );

        // Put beta back
        alpha12L.Set(0,0,beta);
    }
}

template<typename F> 
inline void LUnb( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::LUnb"))
    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Int tHeight = Max(n-1,0);
    t.SetGrid( g );
    t.Resize( tHeight, 1 );

    DistMatrix<F,MC,  STAR> z21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> z1_MR_STAR(g);
    DistMatrix<F,STAR,MC  > a12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > a12_STAR_MR(g);

    for( Int k=0; k<n-1; ++k )
    {
        auto a12      = ViewRange( A, k,   k+1, k+1, n   );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+2 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, n   );
        auto A22      = ViewRange( A, k+1, k+1, n,   n   );
        auto A2       = ViewRange( A, k+1, 0,   n,   n   );

        // Find tau and v such that
        //  |alpha12L a12R| /I - tauP | 1   | | 1 conj(v) |\ = |beta 0|
        //                  \         | v^T |              /
        const F tau = RightReflector( alpha12L, a12R );
        t.Set(k,0,tau);

        // Temporarily set a12 = | 1 v |
        F beta = 0;
        if( alpha12L.IsLocal(0,0) )
            beta = alpha12L.GetLocal(0,0);
        alpha12L.Set(0,0,F(1));

        // A2 := Hous(a12^T,tau)^H A2
        //     = (I - conj(tau) a12^T conj(a12)) A2
        //     = A2 - conj(tau) a12^T (A2^H a12^T)^H
        // -----------------------------------------
        // z1 := A2^H a12^T
        a12_STAR_MC.AlignWith( A2 );
        a12_STAR_MC = a12;
        z1_MR_STAR.AlignWith( A2 );
        Zeros( z1_MR_STAR, n, 1 ); 
        LocalGemv( ADJOINT, F(1), A2, a12_STAR_MC, F(0), z1_MR_STAR );
        z1_MR_STAR.SumOver( A2.ColComm() );
        // A2 := A2 - conj(tau) a12^T z1^H
        LocalGer( -Conj(tau), a12_STAR_MC, z1_MR_STAR, A2 );

        // A22 := A22 Hous(a12^T,tau)
        //      = A22 (I - tau a12^T conj(a12))
        //      = A22 - tau (A22 a12^T) conj(a12)
        // --------------------------------------
        // z21 := A22 a12^T
        a12_STAR_MR.AlignWith( A22 );
        a12_STAR_MR = a12;
        z21_MC_STAR.AlignWith( A22 );
        Zeros( z21_MC_STAR, A22.Height(), 1 );
        LocalGemv( NORMAL, F(1), A22, a12_STAR_MR, F(0), z21_MC_STAR );
        z21_MC_STAR.SumOver( A22.RowComm() );
        // A22 := A22 - tau z21 conj(a12)
        LocalGer( -tau, z21_MC_STAR, a12_STAR_MR, A22 );

        // Put beta back
        if( alpha12L.IsLocal(0,0) )
            alpha12L.SetLocal(0,0,beta);
    }
}

} // namespace hessenberg
} // namespace El

#endif // ifndef EL_HESSENBERG_LUNB_HPP
