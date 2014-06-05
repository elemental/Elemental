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
#ifndef EL_TRMM_LUN_HPP
#define EL_TRMM_LUN_HPP

#include EL_ZEROS_INC

namespace El {
namespace trmm {

template<typename T>
inline void
LocalAccumulateLUN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& U,
  const DistMatrix<T,STAR,MR  >& XTrans,
        DistMatrix<T,MC,  STAR>& Z )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LocalAccumulateLUN");
        if( U.Grid() != XTrans.Grid() ||
            XTrans.Grid() != Z.Grid() )
            LogicError("{U,X,Z} must be distributed over the same grid");
        if( U.Height() != U.Width() ||
            U.Height() != XTrans.Width() ||
            U.Height() != Z.Height() ||
            XTrans.Height() != Z.Width() )
            LogicError
            ("Nonconformal:\n",
             DimsString(U,"U"),"\n",
             DimsString(XTrans,"X'[* ,MR]"),"\n",
             DimsString(Z,"Z[MC,* ]"));
        if( XTrans.RowAlign() != U.RowAlign() ||
            Z.ColAlign() != U.ColAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = Z.Height();
    const Int n = Z.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T> D11(g);

    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto U01 = LockedViewRange( U, 0, k, k,    k+nb );
        auto U11 = LockedViewRange( U, k, k, k+nb, k+nb );

        auto X1Trans = LockedViewRange( XTrans, 0, k, n, k+nb );

        auto Z0 = ViewRange( Z, 0, 0, k,    n );
        auto Z1 = ViewRange( Z, k, 0, k+nb, n );

        D11.AlignWith( U11 );
        D11 = U11;
        MakeTriangular( UPPER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm( NORMAL, orientation, alpha, D11, X1Trans, T(1), Z1 );
        LocalGemm( NORMAL, orientation, alpha, U01, X1Trans, T(1), Z0 );
    }
}

template<typename T>
inline void
LUNA( UnitOrNonUnit diag, const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LUNA");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( U.Height() != U.Width() || U.Width() != X.Height() )
            LogicError
            ("Nonconformal:\n",DimsString(U,"U"),"\n",DimsString(X,"X"));
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T,VR,  STAR> X1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > X1Trans_STAR_MR(g);
    DistMatrix<T,MC,  STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( U );
    X1Trans_STAR_MR.AlignWith( U );
    Z1_MC_STAR.AlignWith( U );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = ViewRange( X, 0, k, m, k+nb );

        X1_VR_STAR = X1;
        X1_VR_STAR.TransposePartialColAllGather( X1Trans_STAR_MR );
        Zeros( Z1_MC_STAR, m, nb );
        LocalAccumulateLUN
        ( TRANSPOSE, diag, T(1), U, X1Trans_STAR_MR, Z1_MC_STAR );

        X1.RowSumScatterFrom( Z1_MC_STAR );
    }
}

template<typename T>
inline void
LUNCOld( UnitOrNonUnit diag, const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LUNCOld");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( U.Height() != U.Width() || U.Width() != X.Height() )
            LogicError
            ("Nonconformal:\n",DimsString(U,"U"),"\n",DimsString(X,"X"));
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        X1_STAR_VR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( LEFT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
 
        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, T(1), X2, U12_STAR_MC, D1Trans_MR_STAR );
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
LUNC( UnitOrNonUnit diag, const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LUNC");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( U.Height() != U.Width() || U.Width() != X.Height() )
            LogicError
            ("Nonconformal:\n",DimsString(U,"U"),"\n",DimsString(X,"X"));
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U01 = LockedViewRange( U, 0, k, k,    k+nb );
        auto U11 = LockedViewRange( U, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01;
        X1Trans_MR_STAR.AlignWith( X0 );
        X1.TransposeColAllGather( X1Trans_MR_STAR );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), U01_MC_STAR, X1Trans_MR_STAR, T(1), X0 );

        U11_STAR_STAR = U11;
        X1_STAR_VR.AlignWith( X1 );
        X1_STAR_VR.TransposePartialRowFilterFrom( X1Trans_MR_STAR );
        LocalTrmm( LEFT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Upper Normal (Non)Unit Trmm
//   X := triu(U)  X, or
//   X := triuu(U) X
template<typename T>
inline void
LUN( UnitOrNonUnit diag, const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("trmm::LUN"))
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Width() )
        LUNA( diag, U, X );
    else
        LUNC( diag, U, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_LUN_HPP
