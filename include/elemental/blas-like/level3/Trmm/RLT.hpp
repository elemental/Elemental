/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRMM_RLT_HPP
#define ELEM_TRMM_RLT_HPP

#include ELEM_AXPY_INC
#include ELEM_MAKETRIANGULAR_INC
#include ELEM_SCALE_INC
#include ELEM_SETDIAGONAL_INC
#include ELEM_TRANSPOSE_INC

#include ELEM_GEMM_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace trmm {

template<typename T>
inline void
LocalAccumulateRLT
( UnitOrNonUnit diag, T alpha,
  const DistMatrix<T>& L,
  const DistMatrix<T,MR,STAR>& XTrans,
        DistMatrix<T,MC,STAR>& ZTrans )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LocalAccumulateRLT");
        if( L.Grid() != XTrans.Grid() ||
            XTrans.Grid() != ZTrans.Grid() )
            LogicError("{L,X,Z} must be distributed over the same grid");
        if( L.Height() != L.Width() ||
            L.Height() != XTrans.Height() ||
            L.Height() != ZTrans.Height() ||
            XTrans.Width() != ZTrans.Width() )
            LogicError
            ("Nonconformal:\n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X^H/T[MR,* ] ~ ",XTrans.Height()," x ",
                                 XTrans.Width(),"\n",
             "  Z^H/T[MC,* ] ~ ",ZTrans.Height()," x ",
                                 ZTrans.Width());
        if( XTrans.ColAlign() != L.RowAlign() ||
            ZTrans.ColAlign() != L.ColAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = ZTrans.Height();
    const Int n = ZTrans.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T> D11(g);

    const Int ratio = Max( g.Height(), g.Width() );

    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1Trans = LockedViewRange( XTrans, k, 0, k+nb, n );

        auto Z1Trans = ViewRange( ZTrans, k,    0, k+nb, n );
        auto Z2Trans = ViewRange( ZTrans, k+nb, 0, m,    n );

        D11.AlignWith( L11 );
        D11 = L11;
        MakeTriangular( LOWER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm( NORMAL, NORMAL, alpha, D11, X1Trans, T(1), Z1Trans );
        LocalGemm( NORMAL, NORMAL, alpha, L21, X1Trans, T(1), Z2Trans );
    }
}

template<typename T>
inline void
RLTA
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::RLTA");
        if( L.Grid() != X.Grid() )
            LogicError("{L,X} must be distributed over the same grid");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);
    DistMatrix<T,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<T,MC,  MR  > Z1Trans(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1Trans_MR_STAR.AlignWith( L );
    Z1Trans_MC_STAR.AlignWith( L );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = ViewRange( X, k, 0, k+nb, n );

        X1.TransposeColAllGather( X1Trans_MR_STAR, conjugate );
        Zeros( Z1Trans_MC_STAR, X1.Width(), X1.Height() );
        LocalAccumulateRLT
        ( diag, T(1), L, X1Trans_MR_STAR, Z1Trans_MC_STAR );

        Z1Trans.RowSumScatterFrom( Z1Trans_MC_STAR );
        Z1Trans_MR_MC.AlignWith( X1 );
        Z1Trans_MR_MC = Z1Trans;
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix(), conjugate );
    }
}

template<typename T>
inline void
RLTC
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::RLTC");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected Adjoint/Transpose option");
        if( L.Height() != L.Width() || X.Width() != L.Height() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrix<T,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);

    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);
         
        auto L10 = LockedViewRange( L, k, 0, k+nb, k    );
        auto L11 = LockedViewRange( L, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, m, k    );
        auto X1 = ViewRange( X, 0, k, m, k+nb );

        X1_VC_STAR = X1;
        L11_STAR_STAR = L11;
        LocalTrmm
        ( RIGHT, LOWER, orientation, diag, T(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        L10Trans_MR_STAR.AlignWith( X0 );
        L10.TransposeColAllGather( L10Trans_MR_STAR, conjugate );
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( NORMAL, NORMAL, T(1), X0, L10Trans_MR_STAR, D1_MC_STAR );
        X1.RowSumScatterUpdate( T(1), D1_MC_STAR );
    }
}

// Right Lower Adjoint/Transpose (Non)Unit Trmm
//   X := X tril(L)^T,
//   X := X tril(L)^H,
//   X := X trilu(L)^T, or
//   X := X trilu(L)^H
template<typename T>
inline void
RLT
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("trmm::RLT"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        RLTA( orientation, diag, L, X );
    else
        RLTC( orientation, diag, L, X );
}

} // namespace trmm
} // namespace elem

#endif // ifndef ELEM_TRMM_RLT_HPP
