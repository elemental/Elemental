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
#ifndef ELEM_TRMM_LUT_HPP
#define ELEM_TRMM_LUT_HPP

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
LocalAccumulateLUT
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T>& U,
  const DistMatrix<T,MC,STAR>& X,
        DistMatrix<T,MR,STAR>& Z )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LocalAccumulateLUT");
        if( U.Grid() != X.Grid() ||
            X.Grid() != Z.Grid() )
            LogicError("{U,X,Z} must be distributed over the same grid");
        if( U.Height() != U.Width() ||
            U.Height() != X.Height() ||
            U.Height() != Z.Height() )
            LogicError
            ("Nonconformal:\n",
             "  U        ~ ",U.Height()," x ",U.Width(),"\n",
             "  X[MC,* ] ~ ",X.Height()," x ",X.Width(),"\n",
             "  Z[MR,* ] ~ ",Z.Height()," x ",Z.Width(),"\n");
        if( X.ColAlign() != U.ColAlign() ||
            Z.ColAlign() != U.RowAlign() )
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

        auto X0 = LockedViewRange( X, 0, 0, k,    n );
        auto X1 = LockedViewRange( X, k, 0, k+nb, n );
 
        auto Z1 = ViewRange( Z, k, 0, k+nb, n );

        D11.AlignWith( U11 );
        D11 = U11;
        MakeTriangular( UPPER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm( orientation, NORMAL, alpha, D11, X1, T(1), Z1 );
        LocalGemm( orientation, NORMAL, alpha, U01, X0, T(1), Z1 );
    }
}

template<typename T>
inline void
LUTA
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LUTA");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( U.Height() != U.Width() || U.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T,MC,STAR> X1_MC_STAR(g);
    DistMatrix<T,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    X1_MC_STAR.AlignWith( U );
    Z1_MR_STAR.AlignWith( U );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = ViewRange( X, 0, k, m, k+nb );

        X1_MC_STAR = X1;
        Zeros( Z1_MR_STAR, m, nb );
        LocalAccumulateLUT
        ( orientation, diag, T(1), U, X1_MC_STAR, Z1_MR_STAR );

        Z1_MR_MC.RowSumScatterFrom( Z1_MR_STAR );
        X1 = Z1_MR_MC;
    }
}

template<typename T>
inline void
LUTCOld
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LUTCOld");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( U.Height() != U.Width() || U.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrix<T,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U01 = LockedViewRange( U, 0, k, k,    k+nb );
        auto U11 = LockedViewRange( U, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        X1_STAR_VR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( LEFT, UPPER, orientation, diag, T(1), U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
        
        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm
        ( orientation, NORMAL, T(1), X0, U01_MC_STAR, D1Trans_MR_STAR );
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
LUTC
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("trmm::LUTC");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( U.Height() != U.Width() || U.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12;
        X1Trans_MR_STAR.AlignWith( X2 ); 
        X1.TransposeColAllGather( X1Trans_MR_STAR );
        LocalGemm
        ( orientation, TRANSPOSE, 
          T(1), U12_STAR_MC, X1Trans_MR_STAR, T(1), X2 );

        U11_STAR_STAR = U11;
        X1_STAR_VR.AlignWith( X1 );
        X1_STAR_VR.TransposePartialRowFilterFrom( X1Trans_MR_STAR );
        LocalTrmm
        ( LEFT, UPPER, orientation, diag, T(1), U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Upper (Conjugate)Transpose (Non)Unit Trmm
//   X := triu(U)^T  X, 
//   X := triu(U)^H  X,
//   X := triuu(U)^T X, or
//   X := triuu(U)^H X
template<typename T>
inline void
LUT
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<T>& U, DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("trmm::LUT"))
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Width() )
        LUTA( orientation, diag, U, X );
    else
        LUTC( orientation, diag, U, X );
}

} // namespace trmm
} // namespace elem

#endif // ifndef ELEM_TRMM_LUT_HPP
