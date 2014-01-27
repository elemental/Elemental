/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_INVERSE_TRIANGULAR_LVAR3_HPP
#define ELEM_INVERSE_TRIANGULAR_LVAR3_HPP

#include ELEM_GEMM_INC
#include ELEM_TRSM_INC

namespace elem {
namespace triang_inv {

template<typename F>
inline void
LVar3Unb( UnitOrNonUnit diag, Matrix<F>& L )
{
    DEBUG_ONLY(
        CallStackEntry cse("triang_inv::LVar3Unb");
        if( L.Height() != L.Width() )
            LogicError("Nonsquare matrices cannot be triangular");
    )
    const Int n = L.Height();
    const Int ldl = L.LDim();
    F* LBuffer = L.Buffer();
    for( Int j=0; j<n; ++j )
    {
        const F lambda = ( diag==NON_UNIT ? LBuffer[j+j*ldl] : F(1) );
        for( Int k=0; k<j; ++k )
            LBuffer[j+k*ldl] /= -lambda;
        blas::Geru
        ( n-(j+1), j, F(1),
          &LBuffer[(j+1)+j*ldl], 1, &LBuffer[j], ldl, 
          &LBuffer[j+1], ldl );
        if( diag == NON_UNIT )
        {
            for( Int k=j+1; k<n; ++k )
                LBuffer[k+j*ldl] /= lambda;
            LBuffer[j+j*ldl] = F(1) / LBuffer[j+j*ldl];
        }
    }
}

template<typename F>
inline void
LVar3( UnitOrNonUnit diag, Matrix<F>& L )
{
    DEBUG_ONLY(
        CallStackEntry cse("triang_inv::LVar3");
        if( L.Height() != L.Width() )
            LogicError("Nonsquare matrices cannot be triangular");
    )
    const Int n = L.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto L10 = ViewRange( L, k,    0, k+nb, k    );
        auto L11 = ViewRange( L, k,    k, k+nb, k+nb );
        auto L20 = ViewRange( L, k+nb, 0, n,    k    );
        auto L21 = ViewRange( L, k+nb, k, n,    k+nb );

        Trsm( LEFT, LOWER, NORMAL, diag, F(-1), L11, L10 );
        Gemm( NORMAL, NORMAL, F(1), L21, L10, F(1), L20 );
        Trsm( RIGHT, LOWER, NORMAL, diag, F(1), L11, L21 );
        LVar3Unb( diag, L11 );
    }
}

template<typename F>
inline void
LVar3( UnitOrNonUnit diag, DistMatrix<F>& L )
{
    DEBUG_ONLY(
        CallStackEntry cse("triang_inv::LVar3");
        if( L.Height() != L.Width() )
            LogicError("Nonsquare matrices cannot be triangular");
    )
    const Grid& g = L.Grid();
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,VC,  STAR> L21_VC_STAR(g);

    const Int n = L.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto L10 = ViewRange( L, k,    0, k+nb, k    );
        auto L11 = ViewRange( L, k,    k, k+nb, k+nb );
        auto L20 = ViewRange( L, k+nb, 0, n,    k    );
        auto L21 = ViewRange( L, k+nb, k, n,    k+nb );

        L10_STAR_VR = L10;
        L11_STAR_STAR = L11;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(-1), L11_STAR_STAR, L10_STAR_VR );

        L21_MC_STAR.AlignWith( L20 );
        L21_MC_STAR = L21;
        L10_STAR_MR.AlignWith( L20 );
        L10_STAR_MR = L10_STAR_VR;
        LocalGemm
        ( NORMAL, NORMAL, F(1), L21_MC_STAR, L10_STAR_MR, F(1), L20 );
        L10 = L10_STAR_MR;

        L21_VC_STAR = L21_MC_STAR;
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, L21_VC_STAR );
        LocalTriangularInverse( LOWER, diag, L11_STAR_STAR );
        L11 = L11_STAR_STAR;
        L21 = L21_VC_STAR;
    }
}

} // namespace triang_inv
} // namespace elem

#endif // ifndef ELEM_INVERSE_TRIANGULAR_LVAR3_HPP
