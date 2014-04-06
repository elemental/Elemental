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
#ifndef ELEM_TRMM_LLT_HPP
#define ELEM_TRMM_LLT_HPP

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
LocalAccumulateLLT
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T>& L,
  const DistMatrix<T,MC,STAR>& X,
        DistMatrix<T,MR,STAR>& Z )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LocalAccumulateLLT");
        if( L.Grid() != X.Grid() || X.Grid() != Z.Grid() )
            LogicError("{L,X,Z} must be distributed over the same grid");
        if( L.Height() != L.Width() ||
            L.Height() != X.Height() ||
            L.Height() != Z.Height() )
            LogicError
            ("Nonconformal:\n",
             "  L        ~ ",L.Height()," x ",L.Width(),"\n",
             "  X[MC,* ] ~ ",X.Height()," x ",X.Width(),"\n",
             "  Z[MR,* ] ~ ",Z.Height()," x ",Z.Width(),"\n");
        if( X.ColAlign() != L.ColAlign() ||
            Z.ColAlign() != L.RowAlign() )
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

        auto X1 = LockedViewRange( X, k,    0, k+nb, n );
        auto X2 = LockedViewRange( X, k+nb, 0, m,    n );

        auto Z1 = ViewRange( Z, k, 0, k+nb, n );

        D11.AlignWith( L11 );
        D11 = L11;
        MakeTriangular( LOWER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm( orientation, NORMAL, alpha, D11, X1, T(1), Z1 );
        LocalGemm( orientation, NORMAL, alpha, L21, X2, T(1), Z1 );
    }
}

template<typename T>
inline void
LLTA
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LLTA");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( L.Height() != L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T,MC,STAR> X1_MC_STAR(g);
    DistMatrix<T,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    X1_MC_STAR.AlignWith( L );
    Z1_MR_STAR.AlignWith( L );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = LockedViewRange( X, 0, k, m, k+nb );

        X1_MC_STAR = X1;
        Zeros( Z1_MR_STAR, m, nb );
        LocalAccumulateLLT
        ( orientation, diag, T(1), L, X1_MC_STAR, Z1_MR_STAR );

        Z1_MR_MC.RowSumScatterFrom( Z1_MR_STAR );
        X1 = Z1_MR_MC;
    }
}
   
template<typename T>
inline void
LLTCOld
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LLTCOld");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( L.Height() != L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        X1_STAR_VR = X1;
        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, orientation, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
 
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm
        ( orientation, NORMAL, T(1), X2, L21_MC_STAR, D1Trans_MR_STAR );
        D1Trans_MR_MC.AlignWith( X1 );
        D1Trans_MR_MC.RowSumScatterFrom( D1Trans_MR_STAR );
        D1.AlignWith( X1 );
        Zeros( D1, nb, n );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix(), conjugate );
        Axpy( T(1), D1, X1 );
    }
}

template<typename T>
inline void
LLTC
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LLTC");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( L.Height() != L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L10 = LockedViewRange( L, k, 0, k+nb, k    );
        auto L11 = LockedViewRange( L, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10;
        X1Trans_MR_STAR.AlignWith( X0 );
        X1.TransposeColAllGather( X1Trans_MR_STAR );
        LocalGemm
        ( orientation, TRANSPOSE, 
          T(1), L10_STAR_MC, X1Trans_MR_STAR, T(1), X0 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.AlignWith( X1 );
        X1_STAR_VR.TransposePartialRowFilterFrom( X1Trans_MR_STAR );
        LocalTrmm
        ( LEFT, LOWER, orientation, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Lower (Conjugate)Transpose (Non)Unit Trmm
//   X := tril(L)^T,
//   X := tril(L)^H,
//   X := trilu(L)^T, or
//   X := trilu(L)^H
template<typename T>
inline void
LLT
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& L, DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("trmm::LLT"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        LLTA( orientation, diag, L, X );
    else
        LLTC( orientation, diag, L, X );
}

} // namespace trmm
} // namespace elem

#endif // ifndef ELEM_TRMM_LLT_HPP
