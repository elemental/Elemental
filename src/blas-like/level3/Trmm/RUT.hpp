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
#ifndef EL_TRMM_RUT_HPP
#define EL_TRMM_RUT_HPP



namespace El {
namespace trmm {

template<typename T>
inline void
LocalAccumulateRUT
( UnitOrNonUnit diag, T alpha,
  const DistMatrix<T>& U,
  const DistMatrix<T,MR,STAR>& XTrans,
        DistMatrix<T,MC,STAR>& ZTrans )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LocalAccumulateRUT");
        if( U.Grid() != XTrans.Grid() ||
            XTrans.Grid() != ZTrans.Grid() )
            LogicError("{U,X,Z} must be distributed over the same grid");
        if( U.Height() != U.Width() ||
            U.Height() != XTrans.Height() ||
            U.Height() != ZTrans.Height() ||
            XTrans.Width() != ZTrans.Width() )
            LogicError
            ("Nonconformal: \n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X^H/T[MR,* ] ~ ",XTrans.Height()," x ",
                                 XTrans.Width(),"\n",
             "  Z^H/T[MC,* ] ~ ",ZTrans.Height()," x ",
                                 ZTrans.Width());
        if( XTrans.ColAlign() != U.RowAlign() ||
            ZTrans.ColAlign() != U.ColAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = ZTrans.Height();
    const Int n = ZTrans.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T> D11(g);

    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto U01 = LockedViewRange( U, 0, k, k,    k+nb );
        auto U11 = LockedViewRange( U, k, k, k+nb, k+nb );

        auto X1Trans = LockedViewRange( XTrans, k, 0, k+nb, n );
   
        auto Z0Trans = ViewRange( ZTrans, 0, 0, k,    n );
        auto Z1Trans = ViewRange( ZTrans, k, 0, k+nb, n );

        D11.AlignWith( U11 );
        D11 = U11;
        MakeTriangular( UPPER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm( NORMAL, NORMAL, alpha, D11, X1Trans, T(1), Z1Trans );
        LocalGemm( NORMAL, NORMAL, alpha, U01, X1Trans, T(1), Z0Trans );
    }
}

template<typename T>
inline void
RUTA
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::RUTA");
        if( U.Grid() != X.Grid() )
            LogicError("{U,X} must be distributed over the same grid");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);
    DistMatrix<T,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<T,MC,  MR  > Z1Trans(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1Trans_MR_STAR.AlignWith( U );
    Z1Trans_MC_STAR.AlignWith( U );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = ViewRange( X, k, 0, k+nb, n );

        X1.TransposeColAllGather( X1Trans_MR_STAR, conjugate );
        Zeros( Z1Trans_MC_STAR, n, nb );
        LocalAccumulateRUT
        ( diag, T(1), U, X1Trans_MR_STAR, Z1Trans_MC_STAR );

        Z1Trans.RowSumScatterFrom( Z1Trans_MC_STAR );
        Z1Trans_MR_MC.AlignWith( X1 );
        Z1Trans_MR_MC = Z1Trans;
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix(), conjugate );
    }
}

template<typename T>
inline void
RUTC
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::RUTC");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected Adjoint/Transpose option");
        if( U.Height() != U.Width() || X.Width() != U.Height() )
            LogicError
            ("Nonconformal: \n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);
    
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, n    );

        auto X1 = ViewRange( X, 0, k,    m, k+nb );
        auto X2 = ViewRange( X, 0, k+nb, m, n    );

        X1_VC_STAR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, orientation, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        U12Trans_MR_STAR.AlignWith( X2 );
        U12.TransposeColAllGather( U12Trans_MR_STAR, conjugate );
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( NORMAL, NORMAL, T(1), X2, U12Trans_MR_STAR, D1_MC_STAR );
        X1.RowSumScatterUpdate( T(1), D1_MC_STAR );
    }
}

// Right Upper Adjoint/Transpose (Non)Unit Trmm
//   X := X triu(U)^T, 
//   X := X triu(U)^H,
//   X := X triuu(U)^T, or
//   X := X triuu(U)^H
template<typename T>
inline void
RUT
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("trmm::RUT"))
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        RUTA( orientation, diag, U, X );
    else
        RUTC( orientation, diag, U, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_RUT_HPP
