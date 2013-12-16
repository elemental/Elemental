/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_TRMM_LUN_HPP
#define ELEM_BLAS_TRMM_LUN_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/blas-like/level1/Transpose.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace internal {

template<typename T>
inline void
LocalTrmmAccumulateLUN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& U,
  const DistMatrix<T,STAR,MR  >& XTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::LocalTrmmAccumulateLUN");
        if( U.Grid() != XTrans_STAR_MR.Grid() ||
            XTrans_STAR_MR.Grid() != Z_MC_STAR.Grid() )
            LogicError("{U,X,Z} must be distributed over the same grid");
        if( U.Height() != U.Width() ||
            U.Height() != XTrans_STAR_MR.Width() ||
            U.Height() != Z_MC_STAR.Height() ||
            XTrans_STAR_MR.Height() != Z_MC_STAR.Width() )
            LogicError
            ("Nonconformal LocalTrmmAccumulateLUN:\n",
             DimsString(U,"U"),"\n",
             DimsString(XTrans_STAR_MR,"X'[* ,MR]"),"\n",
             DimsString(Z_MC_STAR,"Z[MC,* ]"));
        if( XTrans_STAR_MR.RowAlign() != U.RowAlign() ||
            Z_MC_STAR.ColAlign() != U.ColAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MR>
        XLTrans_STAR_MR(g), XRTrans_STAR_MR(g),
        X0Trans_STAR_MR(g), X1Trans_STAR_MR(g), 
        X2Trans_STAR_MR(g);

    DistMatrix<T,MC,STAR>
        ZT_MC_STAR(g),  Z0_MC_STAR(g),
        ZB_MC_STAR(g),  Z1_MC_STAR(g),
                        Z2_MC_STAR(g);

    const Int ratio = Max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    LockedPartitionRight
    ( XTrans_STAR_MR, XLTrans_STAR_MR, XRTrans_STAR_MR, 0 );
    PartitionDown
    ( Z_MC_STAR, ZT_MC_STAR,
                 ZB_MC_STAR, 0 );
    while( UTL.Height() < U.Height() )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        LockedRepartitionRight
        ( XLTrans_STAR_MR, /**/ XRTrans_STAR_MR,
          X0Trans_STAR_MR, /**/ X1Trans_STAR_MR, X2Trans_STAR_MR );

        RepartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
         /**********/ /**********/
                       Z1_MC_STAR,
          ZB_MC_STAR,  Z2_MC_STAR );

        D11.AlignWith( U11 );
        //--------------------------------------------------------------------//
        D11 = U11;
        MakeTriangular( UPPER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm
        ( NORMAL, orientation, alpha, D11, X1Trans_STAR_MR, T(1), Z1_MC_STAR );
        LocalGemm
        ( NORMAL, orientation, alpha, U01, X1Trans_STAR_MR, T(1), Z0_MC_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlideLockedPartitionRight
        ( XLTrans_STAR_MR,                  /**/ XRTrans_STAR_MR,
          X0Trans_STAR_MR, X1Trans_STAR_MR, /**/ X2Trans_STAR_MR );

        SlidePartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
                       Z1_MC_STAR,
         /**********/ /**********/
          ZB_MC_STAR,  Z2_MC_STAR );
    }
    PopBlocksizeStack();
}

template<typename T>
inline void
TrmmLUNA
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmULNA");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( U.Height() != U.Width() || U.Width() != X.Height() )
            LogicError
            ("Nonconformal TrmmLUNA:\n",
             DimsString(U,"U"),"\n",DimsString(X,"X"));
    )
    const Grid& g = U.Grid();

    DistMatrix<T>
        XL(g), XR(g),
        X0(g), X1(g), X2(g);

    DistMatrix<T,VR,  STAR> X1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > X1Trans_STAR_MR(g);
    DistMatrix<T,MC,  STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( U );
    X1Trans_STAR_MR.AlignWith( U );
    Z1_MC_STAR.AlignWith( U );

    PartitionRight( X, XL, XR, 0 );
    while( XL.Width() < X.Width() )
    {
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        //--------------------------------------------------------------------//
        X1_VR_STAR = X1;
        X1Trans_STAR_MR.TransposeFrom( X1_VR_STAR );
        Zeros( Z1_MC_STAR, X1.Height(), X1.Width() );
        LocalTrmmAccumulateLUN
        ( TRANSPOSE, diag, alpha, U, X1Trans_STAR_MR, Z1_MC_STAR );

        X1.SumScatterFrom( Z1_MC_STAR );
        //--------------------------------------------------------------------//

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
}

template<typename T>
inline void
TrmmLUNCOld
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmLUNCOld");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( U.Height() != U.Width() || U.Width() != X.Height() )
            LogicError
            ("Nonconformal TrmmLUN:\n",
             DimsString(U,"U"),"\n",DimsString(X,"X"));
    )
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);
    DistMatrix<T> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,   U00, /**/ U01, U02,
         /*************/  /******************/
               /**/        U10, /**/ U11, U12,
          UBL, /**/ UBR,   U20, /**/ U21, U22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        U12_STAR_MC.AlignWith( X2 );
        D1Trans_MR_STAR.AlignWith( X1 );
        D1Trans_MR_MC.AlignWith( X1 );
        D1.AlignWith( X1 );
        //--------------------------------------------------------------------//
        X1_STAR_VR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( LEFT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
 
        U12_STAR_MC = U12;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, T(1), X2, U12_STAR_MC, D1Trans_MR_STAR );
        D1Trans_MR_MC.SumScatterFrom( D1Trans_MR_STAR );
        Zeros( D1, X1.Height(), X1.Width() );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix() );
        Axpy( T(1), D1, X1 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
}

template<typename T>
inline void
TrmmLUNC
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmLUNC");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( U.Height() != U.Width() || U.Width() != X.Height() )
            LogicError
            ("Nonconformal TrmmLUN:\n",
             DimsString(U,"U"),"\n",DimsString(X,"X"));
    )
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);
    DistMatrix<T> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,   U00, /**/ U01, U02,
         /*************/  /******************/
               /**/        U10, /**/ U11, U12,
          UBL, /**/ UBR,   U20, /**/ U21, U22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        U01_MC_STAR.AlignWith( X0 );
        X1Trans_MR_STAR.AlignWith( X0 );
        X1_STAR_VR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        U01_MC_STAR = U01;
        X1Trans_MR_STAR.TransposeFrom( X1 );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), U01_MC_STAR, X1Trans_MR_STAR, T(1), X0 );

        U11_STAR_STAR = U11;
        X1_STAR_VR.TransposeFrom( X1Trans_MR_STAR );
        LocalTrmm( LEFT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
}

// Left Upper Normal (Non)Unit Trmm
//   X := triu(U)  X, or
//   X := triuu(U) X
template<typename T>
inline void
TrmmLUN
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("internal::TrmmLUN"))
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Width() )
        TrmmLUNA( diag, alpha, U, X );
    else
        TrmmLUNC( diag, alpha, U, X );
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_BLAS_TRMM_LUN_HPP
