/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CHOLESKY_LVAR2_HPP
#define EL_CHOLESKY_LVAR2_HPP

#include "./LVar3.hpp"

// TODO: Reverse variants

namespace El {
namespace cholesky {

template<typename F> 
inline void
LVar2( Matrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::LVar2");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto A10 = ViewRange( A, k,    0, k+nb, k    );
        auto A11 = ViewRange( A, k,    k, k+nb, k+nb );
        auto A20 = ViewRange( A, k+nb, 0, n,    k    );
        auto A21 = ViewRange( A, k+nb, k, n,    k+nb );

        Herk( LOWER, NORMAL, F(-1), A10, F(1), A11 );
        cholesky::LVar3Unb( A11 );
        Gemm( NORMAL, ADJOINT, F(-1), A20, A10, F(1), A21 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11, A21 );
    }
}

template<typename F> 
inline void
LVar2( DistMatrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::LVar2");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& g = A.Grid();
    DistMatrix<F,MR,  STAR> A10Adj_MR_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,MC,  STAR> X11_MC_STAR(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto A10 = ViewRange( A, k,    0, k+nb, k    );
        auto A11 = ViewRange( A, k,    k, k+nb, k+nb );
        auto A20 = ViewRange( A, k+nb, 0, n,    k    );
        auto A21 = ViewRange( A, k+nb, k, n,    k+nb );
 
        A10Adj_MR_STAR.AlignWith( A10 );
        A10.AdjointColAllGather( A10Adj_MR_STAR );
        X11_MC_STAR.AlignWith( A10 );
        LocalGemm( NORMAL, NORMAL, F(1), A10, A10Adj_MR_STAR, X11_MC_STAR );
        A11.RowSumScatterUpdate( F(-1), X11_MC_STAR );

        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        X21_MC_STAR.AlignWith( A20 );
        LocalGemm( NORMAL, NORMAL, F(1), A20, A10Adj_MR_STAR, X21_MC_STAR );
        A21.RowSumScatterUpdate( F(-1), X21_MC_STAR );

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_LVAR2_HPP
