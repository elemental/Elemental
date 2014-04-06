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
#ifndef ELEM_TRMM_LLN_HPP
#define ELEM_TRMM_LLN_HPP

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
LocalAccumulateLLN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& L,
  const DistMatrix<T,STAR,MR  >& XTrans,
        DistMatrix<T,MC,  STAR>& Z )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LocalAccumulateLLN");
        if( L.Grid() != XTrans.Grid() ||
            XTrans.Grid() != Z.Grid() )
            LogicError("{L,X,Z} must be distributed over the same grid");
        if( L.Height() != L.Width() ||
            L.Height() != XTrans.Width() ||
            L.Height() != Z.Height() ||
            XTrans.Height() != Z.Width() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X^H/T[* ,MR] ~ ",XTrans.Height()," x ",
                                 XTrans.Width(),"\n",
             "  Z[MC,* ] ~ ",Z.Height()," x ",Z.Width());
        if( XTrans.RowAlign() != L.RowAlign() ||
            Z.ColAlign() != L.ColAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = Z.Height();
    const Int n = Z.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T> D11(g);
    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1Trans = LockedViewRange( XTrans, 0, k, n, k+nb );

        auto Z1 = ViewRange( Z, k,    0, k+nb, n );
        auto Z2 = ViewRange( Z, k+nb, 0, m,    n );

        D11.AlignWith( L11 );
        D11 = L11;
        MakeTriangular( LOWER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm( NORMAL, orientation, alpha, D11, X1Trans, T(1), Z1 );
        LocalGemm( NORMAL, orientation, alpha, L21, X1Trans, T(1), Z2 );
    }
}

template<typename T>
inline void
LLNA( UnitOrNonUnit diag, const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LLNA");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() != L.Width() || L.Width() != X.Height() )
            LogicError
            ("Nonconformal: \n"
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T,VR,  STAR> X1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > X1Trans_STAR_MR(g);
    DistMatrix<T,MC,  STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( L );
    X1Trans_STAR_MR.AlignWith( L );
    Z1_MC_STAR.AlignWith( L );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = ViewRange( X, 0, k, m, k+nb );

        X1_VR_STAR = X1;
        X1_VR_STAR.TransposePartialColAllGather( X1Trans_STAR_MR );
        Zeros( Z1_MC_STAR, m, nb );
        LocalAccumulateLLN
        ( TRANSPOSE, diag, T(1), L, X1Trans_STAR_MR, Z1_MC_STAR );
        X1.RowSumScatterFrom( Z1_MC_STAR );
    }
}

template<typename T>
inline void
LLNCOld( UnitOrNonUnit diag, const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LLNCOld");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() != L.Width() || L.Width() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L10 = LockedViewRange( L, k, 0, k+nb, k    );
        auto L11 = LockedViewRange( L, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        L11_STAR_STAR = L11;
        X1_STAR_VR = X1;
        LocalTrmm( LEFT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;

        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, T(1), X0, L10_STAR_MC, D1Trans_MR_STAR );
        D1Trans_MR_MC.AlignWith( X1 );
        D1Trans_MR_MC.RowSumScatterFrom( D1Trans_MR_STAR );
        D1.AlignWith( X1 );
        Zeros( D1, nb, n );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix() );
        Axpy( T(1), D1, X1 );
    }
}

template<typename T>
inline void
LLNC( UnitOrNonUnit diag, const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LLNC");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() != L.Width() || L.Width() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;
        X1Trans_MR_STAR.AlignWith( X2 );
        X1.TransposeColAllGather( X1Trans_MR_STAR );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), L21_MC_STAR, X1Trans_MR_STAR, T(1), X2 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.AlignWith( X1 );
        X1_STAR_VR.TransposePartialRowFilterFrom( X1Trans_MR_STAR );
        LocalTrmm( LEFT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Lower Normal (Non)Unit Trmm 
//   X := tril(L)  X, or
//   X := trilu(L) X
template<typename T>
inline void
LLN( UnitOrNonUnit diag, const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("trmm::LLN"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        LLNA( diag, L, X );
    else
        LLNC( diag, L, X );
}

} // namespace trmm
} // namespace elem

#endif // ifndef ELEM_TRMM_LLN_HPP
