/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HESSENBERG_U_HPP
#define ELEM_HESSENBERG_U_HPP

#include "./UUnb.hpp"
#include "./PanelU.hpp"

namespace elem {
namespace hessenberg {

template<typename F>
inline void U( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("hessenberg::U");
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
        hessenberg::PanelU( ABR, t1, UB1, VB1, G11 );

        auto A0R = ViewRange( A,   0,  k,    k,   n  );
        auto AB2 = ViewRange( A,   k,  k+nb, n,   n  );
        auto U21 = ViewRange( UB1, nb, 0,    n-k, nb );
        auto V21 = ViewRange( VB1, nb, 0,    n-k, nb );

        // A0R := A0R - ((A0R UB1) inv(G11)^H) UB1^H
        // -----------------------------------------
        Gemm( NORMAL, NORMAL, F(1), A0R, UB1, V01 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), G11, V01 );
        Gemm( NORMAL, ADJOINT, F(-1), V01, UB1, F(1), A0R );
            
        // AB2 := (I - UB1 inv(G11) UB1^H)(AB2 - VB1 inv(G11)^H U21^H)
        // -----------------------------------------------------------
        // AB2 := AB2 - VB1 inv(G11)^H U21^H 
        // (note: VB1 is overwritten) 
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), G11, VB1 );
        Gemm( NORMAL, ADJOINT, F(-1), VB1, U21, F(1), AB2 );
        // AB2 := AB2 - UB1 (inv(G11) (UB1^H AB2))
        //      = AB2 - UB1 ((AB2^H UB1) inv(G11)^H)^H
        // (note: V21 is used as scratch space)
        Gemm( ADJOINT, NORMAL, F(1), AB2, UB1, F(0), V21 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), G11, V21 );
        Gemm( NORMAL, ADJOINT, F(-1), UB1, V21, F(1), AB2 );
    }
}

template<typename F> 
inline void U( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::U"))
    // TODO: Blocked algorithm
    UUnb( A, t );
}

} // namespace hessenberg
} // namespace elem

#endif // ifndef ELEM_HESSENBERG_U_HPP
