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
#ifndef ELEM_TRMM_RUN_HPP
#define ELEM_TRMM_RUN_HPP

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
LocalAccumulateRUN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& U,
  const DistMatrix<T,STAR,MC  >& X,
        DistMatrix<T,MR,  STAR>& ZTrans )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LocalAccumulateRUN");
        if( U.Grid() != X.Grid() ||
            X.Grid() != ZTrans.Grid() )
            LogicError("{U,X,Z} must be distributed over the same grid");
        if( U.Height() != U.Width() ||
            U.Height() != X.Width() ||
            U.Height() != ZTrans.Height() )
            LogicError
            ("Nonconformal:\n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X[* ,MC] ~ ",X.Height()," x ",X.Width(),"\n",
             "  Z^H/T[MR,* ] ~ ",ZTrans.Height()," x ",
                                 ZTrans.Width());
        if( X.RowAlign() != U.ColAlign() ||
            ZTrans.ColAlign() != U.RowAlign() )
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

        auto X0 = LockedViewRange( X, 0, 0, n, k    );
        auto X1 = LockedViewRange( X, 0, k, n, k+nb );
        
        auto Z1Trans = ViewRange( ZTrans, k, 0, k+nb, n );

        D11.AlignWith( U11 );
        D11 = U11;
        MakeTriangular( UPPER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm( orientation, orientation, alpha, D11, X1, T(1), Z1Trans );
        LocalGemm( orientation, orientation, alpha, U01, X0, T(1), Z1Trans );
    }
}

template<typename T>
inline void
RUNA( UnitOrNonUnit diag, const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::RUNA");
        if( U.Grid() != X.Grid() )
            LogicError("{U,X} must be distributed over the same grid");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T,STAR,VC  > X1_STAR_VC(g);
    DistMatrix<T,STAR,MC  > X1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1_STAR_VC.AlignWith( U );
    X1_STAR_MC.AlignWith( U );
    Z1Trans_MR_STAR.AlignWith( U );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = ViewRange( X, k, 0, k+nb, n );

        X1_STAR_VC = X1;
        X1_STAR_MC = X1_STAR_VC;
        Zeros( Z1Trans_MR_STAR, X1.Width(), X1.Height() );
        LocalAccumulateRUN
        ( TRANSPOSE, diag, T(1), U, X1_STAR_MC, Z1Trans_MR_STAR );

        Z1Trans_MR_MC.AlignWith( X1 );
        Z1Trans_MR_MC.RowSumScatterFrom( Z1Trans_MR_STAR );
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix() );
    }
}

template<typename T>
inline void
RUNCOld( UnitOrNonUnit diag, const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::RUNCOld");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( U.Height() != U.Width() || X.Width() != U.Height() )
            LogicError
            ("Nonconformal:\n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T,MR,  STAR> U01_MR_STAR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);    
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);
    
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U01 = LockedViewRange( U, 0, k, k,    k+nb );
        auto U11 = LockedViewRange( U, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, m, k    );
        auto X1 = ViewRange( X, 0, k, m, k+nb );

        X1_VC_STAR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        U01_MR_STAR.AlignWith( X0 );
        U01_MR_STAR = U01;
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( NORMAL, NORMAL, T(1), X0, U01_MR_STAR, D1_MC_STAR );
        X1.RowSumScatterUpdate( T(1), D1_MC_STAR );
    }
}

template<typename T>
inline void
RUNC( UnitOrNonUnit diag, const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::RUNC");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( U.Height() != U.Width() || X.Width() != U.Height() )
            LogicError
            ("Nonconformal:\n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> X1_MC_STAR(g);
    
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, n    );

        auto X1 = ViewRange( X, 0, k,    m, k+nb );
        auto X2 = ViewRange( X, 0, k+nb, m, n    );

        X1_MC_STAR.AlignWith( X2 );
        X1_MC_STAR = X1;
        U12Trans_MR_STAR.AlignWith( X2 );
        U12.TransposeColAllGather( U12Trans_MR_STAR );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), X1_MC_STAR, U12Trans_MR_STAR, T(1), X2 );

        U11_STAR_STAR = U11;
        X1_VC_STAR.AlignWith( X1 );
        X1_VC_STAR = X1_MC_STAR;
        LocalTrmm
        ( RIGHT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
    }
}

// Right Upper Normal (Non)Unit Trmm
//   X := X triu(U), and
//   X := X triuu(U)
template<typename T>
inline void
RUN( UnitOrNonUnit diag, const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("trmm::RUN"))
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        RUNA( diag, U, X );
    else
        RUNC( diag, U, X );
}

} // namespace trmm
} // namespace elem

#endif // ifndef ELEM_TRMM_RUN_HPP
