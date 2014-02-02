/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HESSENBERG_L_HPP
#define ELEM_HESSENBERG_L_HPP

#include "./LUnb.hpp"
#include "./PanelL.hpp"

namespace elem {
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
        hessenberg::PanelL( ABR, t1, UB1, VB1, G11 );

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
    DEBUG_ONLY(CallStackEntry cse("hessenberg::L"))
    // TODO: Blocked algorithm
    LUnb( A, t );
}

} // namespace hessenberg
} // namespace elem

#endif // ifndef ELEM_HESSENBERG_L_HPP
