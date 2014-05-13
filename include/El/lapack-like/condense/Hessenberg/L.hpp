/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HESSENBERG_L_HPP
#define EL_HESSENBERG_L_HPP

#include "./LUnb.hpp"
#include "./LPan.hpp"

namespace El {
namespace hessenberg {

template<typename F>
inline void L( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("hessenberg::L");
        // Is this requirement necessary?!?
        if( t.Viewing() )
            LogicError("t must not be a view");
    )
    const Int n = A.Height();
    t.Resize( Max(n-1,0), 1 );

    Matrix<F> UB1, V01, VB1, G11;

    const Int bsize = Blocksize();
    for( Int k=0; k<n-1; k+=bsize )
    {
        const Int nb = Min(bsize,n-1-k);
        auto ABR = ViewRange( A, k,    k,    n, n );
        auto A22 = ViewRange( A, k+nb, k+nb, n, n );

        auto t1 = View( t, k, 0, nb, 1 );
        UB1.Resize( n-k, nb );
        VB1.Resize( n-k, nb );
        G11.Resize( nb,  nb );
        hessenberg::LPan( ABR, t1, UB1, VB1, G11 );

        auto AB0 = ViewRange( A,   k,    0, n,   k  );
        auto A2R = ViewRange( A,   k+nb, k, n,   n  );
        auto U21 = ViewRange( UB1, nb,   0, n-k, nb );
        auto V21 = ViewRange( VB1, nb,   0, n-k, nb );

        // AB0 := AB0 - (UB1 inv(G11)^H UB1^H AB0)
        //      = AB0 - (UB1 ((AB0^H UB1) inv(G11))^H)
        // -------------------------------------------
        Gemm( ADJOINT, NORMAL, F(1), AB0, UB1, V01 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11, V01 ); 
        Gemm( NORMAL, ADJOINT, F(-1), UB1, V01, F(1), AB0 );

        // A2R := (A2R - U21 inv(G11)^H VB1^H)(I - UB1 inv(G11) UB1^H)
        // -----------------------------------------------------------
        // A2R := A2R - U21 inv(G11)^H VB1^H
        // (note: VB1 is overwritten)
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11, VB1 );
        Gemm( NORMAL, ADJOINT, F(-1), U21, VB1, F(1), A2R );
        // A2R := A2R - ((A2R UB1) inv(G11)) UB1^H
        Gemm( NORMAL, NORMAL, F(1), A2R, UB1, F(0), V21 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11, V21 );
        Gemm( NORMAL, ADJOINT, F(-1), V21, UB1, F(1), A2R );
    }
}

template<typename F> 
inline void L( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("hessenberg::L");
        if( A.Grid() != t.Grid() )
            LogicError("A and t must be distributed over the same grid");
        // Is this requirement necessary?!?
        if( t.Viewing() )
            LogicError("t must not be a view");
    )
    const Grid& g = A.Grid();
    const Int n = A.Height();
    t.Resize( Max(n-1,0), 1 );

    DistMatrix<F,MC,STAR> UB1_MC_STAR(g), V21_MC_STAR(g);
    DistMatrix<F,MR,STAR> V01_MR_STAR(g), VB1_MR_STAR(g), UB1_MR_STAR(g);
    DistMatrix<F,STAR,STAR> G11_STAR_STAR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n-1; k+=bsize )
    {
        const Int nb = Min(bsize,n-1-k);
        auto ABR = ViewRange( A, k,    k,    n, n );
        auto A22 = ViewRange( A, k+nb, k+nb, n, n );

        auto t1 = View( t, k, 0, nb, 1 );
        UB1_MC_STAR.AlignWith( ABR );
        UB1_MR_STAR.AlignWith( ABR );
        VB1_MR_STAR.AlignWith( ABR );
        UB1_MC_STAR.Resize( n-k, nb );
        UB1_MR_STAR.Resize( n-k, nb );
        VB1_MR_STAR.Resize( n-k, nb );
        G11_STAR_STAR.Resize( nb, nb );
        hessenberg::LPan
        ( ABR, t1, UB1_MC_STAR, UB1_MR_STAR, VB1_MR_STAR, G11_STAR_STAR );

        auto AB0 = ViewRange( A,   k,    0, n,   k  );
        auto A2R = ViewRange( A,   k+nb, k, n,   n  );

        auto U21_MC_STAR = LockedViewRange( UB1_MC_STAR, nb, 0, n-k, nb );

        // AB0 := AB0 - (UB1 inv(G11)^H UB1^H AB0)
        //      = AB0 - (UB1 ((AB0^H UB1) inv(G11))^H)
        // -------------------------------------------
        V01_MR_STAR.AlignWith( AB0 );
        Zeros( V01_MR_STAR, k, nb );
        LocalGemm( ADJOINT, NORMAL, F(1), AB0, UB1_MC_STAR, F(0), V01_MR_STAR );
        V01_MR_STAR.SumOver( AB0.ColComm() );
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11_STAR_STAR, V01_MR_STAR );
        LocalGemm
        ( NORMAL, ADJOINT, F(-1), UB1_MC_STAR, V01_MR_STAR, F(1), AB0 );

        // A2R := (A2R - U21 inv(G11)^H VB1^H)(I - UB1 inv(G11) UB1^H)
        // -----------------------------------------------------------
        // A2R := A2R - U21 inv(G11)^H VB1^H
        // (note: VB1 is overwritten)
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11_STAR_STAR, VB1_MR_STAR );
        LocalGemm
        ( NORMAL, ADJOINT, F(-1), U21_MC_STAR, VB1_MR_STAR, F(1), A2R );
        // A2R := A2R - ((A2R UB1) inv(G11)) UB1^H
        V21_MC_STAR.AlignWith( A2R );
        Zeros( V21_MC_STAR, A2R.Height(), nb );
        LocalGemm( NORMAL, NORMAL, F(1), A2R, UB1_MR_STAR, F(0), V21_MC_STAR );
        V21_MC_STAR.SumOver( A2R.RowComm() );
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11_STAR_STAR, V21_MC_STAR );
        LocalGemm
        ( NORMAL, ADJOINT, F(-1), V21_MC_STAR, UB1_MR_STAR, F(1), A2R );
    }
}

} // namespace hessenberg
} // namespace El

#endif // ifndef EL_HESSENBERG_L_HPP
