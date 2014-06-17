/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HESSENBERG_UUNB_HPP
#define EL_HESSENBERG_UUNB_HPP

namespace El {
namespace hessenberg {

template<typename F>
inline void UUnb( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::UUnb"))
    const Int n = A.Height();
    const Int tHeight = Max(n-1,0);
    t.Resize( tHeight, 1 );

    // Temporary products
    Matrix<F> x1, x12Adj;

    for( Int k=0; k<n-1; ++k )
    {
        auto a21      = ViewRange( A, k+1, k,   n,   k+1 );
        auto alpha21T = ViewRange( A, k+1, k,   k+2, k+1 );
        auto a21B     = ViewRange( A, k+2, k,   n,   k+1 );
        auto A22      = ViewRange( A, k+1, k+1, n,   n   );
        auto A2       = ViewRange( A, 0,   k+1, n,   n   );

        // Find tau and v such that
        //  / I - tau | 1 | | 1, v^H | \ | alpha21T | = | beta |
        //  \         | v |            / |     a21B |   |    0 |
        const F tau = LeftReflector( alpha21T, a21B );
        t.Set(k,0,tau);

        // Temporarily set a21 := | 1 |
        //                        | v |
        const F beta = alpha21T.Get(0,0);
        alpha21T.Set(0,0,F(1));

        // A2 := A2 Hous(a21,tau)^H
        //     = A2 (I - conj(tau) a21 a21^H)
        //     = A2 - conj(tau) (A2 a21) a21^H
        // -----------------------------------
        // x1 := A2 a21
        Zeros( x1, n, 1 );
        Gemv( NORMAL, F(1), A2, a21, F(0), x1 );
        // A2 := A2 - conj(tau) x1 a21^H
        Ger( -Conj(tau), x1, a21, A2 ); 

        // A22 := Hous(a21,tau) A22
        //      = (I - tau a21 a21^H) A22
        //      = A22 - tau a21 (A22^H a21)^H
        // ----------------------------------
        // x12^H := (a21^H A22)^H = A22^H a21
        Zeros( x12Adj, A22.Width(), 1 );
        Gemv( ADJOINT, F(1), A22, a21, F(0), x12Adj );
        // A22 := A22 - tau a21 x12
        Ger( -tau, a21, x12Adj, A22 );

        // Put beta back
        alpha21T.Set(0,0,beta);
    }
}

template<typename F> 
inline void UUnb( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::UUnb"))
    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Int tHeight = Max(n-1,0);
    t.SetGrid( g );
    t.Resize( tHeight, 1 );

    DistMatrix<F,MC,STAR> a21_MC_STAR(g);
    DistMatrix<F,MR,STAR> a21_MR_STAR(g);
    DistMatrix<F,MC,STAR> x1_MC_STAR(g); 
    DistMatrix<F,MR,STAR> x12Adj_MR_STAR(g);

    for( Int k=0; k<n-1; ++k )
    {
        auto a21      = ViewRange( A, k+1, k,   n,   k+1 );
        auto alpha21T = ViewRange( A, k+1, k,   k+2, k+1 );
        auto a21B     = ViewRange( A, k+2, k,   n,   k+1 );
        auto A22      = ViewRange( A, k+1, k+1, n,   n   );
        auto A2       = ViewRange( A, 0,   k+1, n,   n   );

        // Find tau and v such that
        //  / I - tau | 1 | | 1, v^H | \ | alpha21T | = | beta |
        //  \         | v |            / |     a21B |   |    0 |
        const F tau = LeftReflector( alpha21T, a21B );
        t.Set(k,0,tau);

        // Temporarily set a21 := | 1 |
        //                        | v |
        F beta=0;
        if( alpha21T.IsLocal(0,0) )
            beta = alpha21T.GetLocal(0,0);
        alpha21T.Set(0,0,F(1));

        // A2 := A2 Hous(a21,tau)^H
        //     = A2 (I - conj(tau) a21 a21^H)
        //     = A2 - conj(tau) (A2 a21) a21^H
        // -----------------------------------
        // x1 := A2 a21
        a21_MR_STAR.AlignWith( A2 );
        a21_MR_STAR = a21;
        x1_MC_STAR.AlignWith( A2 );
        Zeros( x1_MC_STAR, n, 1 );
        LocalGemv( NORMAL, F(1), A2, a21_MR_STAR, F(0), x1_MC_STAR );
        x1_MC_STAR.SumOver( A2.RowComm() );
        // A2 := A2 - conj(tau) x1 a21^H
        LocalGer( -Conj(tau), x1_MC_STAR, a21_MR_STAR, A2 ); 

        // A22 := Hous(a21,tau) A22
        //      = (I - tau a21 a21^H) A22
        //      = A22 - tau a21 (A22^H a21)^H
        // ----------------------------------
        // x12^H := (a21^H A22)^H = A22^H a21
        a21_MC_STAR.AlignWith( A22 );
        a21_MC_STAR = a21;
        x12Adj_MR_STAR.AlignWith( A22 );
        Zeros( x12Adj_MR_STAR, A22.Width(), 1 );
        LocalGemv( ADJOINT, F(1), A22, a21_MC_STAR, F(0), x12Adj_MR_STAR );
        x12Adj_MR_STAR.SumOver( A22.ColComm() );
        // A22 := A22 - tau a21 x12
        LocalGer( -tau, a21_MC_STAR, x12Adj_MR_STAR, A22 );

        // Put beta back 
        if( alpha21T.IsLocal(0,0) )
            alpha21T.SetLocal(0,0,beta);
    }
}

} // namespace hessenberg
} // namespace El

#endif // ifndef EL_HESSENBERG_UUNB_HPP
