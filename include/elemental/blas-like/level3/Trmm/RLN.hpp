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
#ifndef ELEM_BLAS_TRMM_RLN_HPP
#define ELEM_BLAS_TRMM_RLN_HPP

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
LocalTrmmAccumulateRLN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& L,
  const DistMatrix<T,STAR,MC  >& X_STAR_MC,
        DistMatrix<T,MR,  STAR>& ZTrans_MR_STAR )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::LocalTrmmAccumulateRLN");
        if( L.Grid() != X_STAR_MC.Grid() || 
            X_STAR_MC.Grid() != ZTrans_MR_STAR.Grid() )
            LogicError("{L,X,Z} must be distributed over the same grid");
        if( L.Height() != L.Width() ||
            L.Height() != X_STAR_MC.Width() ||
            L.Height() != ZTrans_MR_STAR.Height() )
            LogicError
            ("Nonconformal LocalTrmmAccumulateRLN:\n",
             DimsString(L,"L"),"\n",
             DimsString(X_STAR_MC,"X[* ,MC]"),"\n",
             DimsString(ZTrans_MR_STAR,"Z'[MR,* ]"));
        if( X_STAR_MC.RowAlign() != L.ColAlign() ||
            ZTrans_MR_STAR.ColAlign() != L.RowAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MC>
        XL_STAR_MC(g), XR_STAR_MC(g),
        X0_STAR_MC(g), X1_STAR_MC(g), X2_STAR_MC(g);

    DistMatrix<T,MR,STAR>
        ZTTrans_MR_STAR(g),  Z0Trans_MR_STAR(g),
        ZBTrans_MR_STAR(g),  Z1Trans_MR_STAR(g),
                             Z2Trans_MR_STAR(g);

    const Int ratio = Max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionRight( X_STAR_MC,  XL_STAR_MC, XR_STAR_MC, 0 );
    PartitionDown
    ( ZTrans_MR_STAR, ZTTrans_MR_STAR,
                      ZBTrans_MR_STAR, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionRight
        ( XL_STAR_MC, /**/ XR_STAR_MC,
          X0_STAR_MC, /**/ X1_STAR_MC, X2_STAR_MC );

        RepartitionDown
        ( ZTTrans_MR_STAR,  Z0Trans_MR_STAR,
         /***************/ /***************/
                            Z1Trans_MR_STAR,
          ZBTrans_MR_STAR,  Z2Trans_MR_STAR );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        MakeTriangular( LOWER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm
        ( orientation, orientation,
          alpha, D11, X1_STAR_MC, T(1), Z1Trans_MR_STAR );
        LocalGemm
        ( orientation, orientation,
          alpha, L21, X2_STAR_MC, T(1), Z1Trans_MR_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlideLockedPartitionRight
        ( XL_STAR_MC,             /**/ XR_STAR_MC,
          X0_STAR_MC, X1_STAR_MC, /**/ X2_STAR_MC );

        SlidePartitionDown
        ( ZTTrans_MR_STAR,  Z0Trans_MR_STAR,
                            Z1Trans_MR_STAR,
         /***************/ /***************/
          ZBTrans_MR_STAR,  Z2Trans_MR_STAR );
    }
    PopBlocksizeStack();
}

template<typename T>
inline void
TrmmRLNA
( UnitOrNonUnit diag, 
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmRLNA");
        if( L.Grid() != X.Grid() )
            LogicError("{L,X} must be distributed over the same grid");
    )
    const Grid& g = L.Grid();

    DistMatrix<T>
        XT(g),  X0(g),
        XB(g),  X1(g),
                X2(g);

    DistMatrix<T,STAR,VC  > X1_STAR_VC(g);
    DistMatrix<T,STAR,MC  > X1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1_STAR_VC.AlignWith( L );
    X1_STAR_MC.AlignWith( L );
    Z1Trans_MR_STAR.AlignWith( L );

    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XT.Height() < X.Height() )
    {
        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        Z1Trans_MR_MC.AlignWith( X1 );
        //--------------------------------------------------------------------//
        X1_STAR_VC = X1;
        X1_STAR_MC = X1_STAR_VC;
        Zeros( Z1Trans_MR_STAR, X1.Width(), X1.Height() );
        LocalTrmmAccumulateRLN
        ( TRANSPOSE, diag, alpha, L, X1_STAR_MC, Z1Trans_MR_STAR );

        Z1Trans_MR_MC.SumScatterFrom( Z1Trans_MR_STAR );
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix() );
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
}
    
template<typename T>
inline void
TrmmRLNCOld
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmRLNCOld");
        if( L.Grid() != X.Grid() )
            LogicError
            ("L and X must be distributed over the same grid");
        if( L.Height() != L.Width() || X.Width() != L.Height() )
            LogicError
            ("Nonconformal TrmmRLNC:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> L21_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );
 
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        L21_MR_STAR.AlignWith( X2 );
        D1_MC_STAR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        X1_VC_STAR = X1;
        L11_STAR_STAR = L11;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        L21_MR_STAR = L21;
        LocalGemm( NORMAL, NORMAL, T(1), X2, L21_MR_STAR, D1_MC_STAR );
        X1.SumScatterUpdate( T(1), D1_MC_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
}

template<typename T>
inline void
TrmmRLNC
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmRLNC");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() != L.Width() || X.Width() != L.Height() )
            LogicError
            ("Nonconformal TrmmRLNC:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
    )
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> X1_MC_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );
 
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        X1_MC_STAR.AlignWith( X0 );
        L10Trans_MR_STAR.AlignWith( X0 );
        X1_VC_STAR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        X1_MC_STAR = X1;
        L10Trans_MR_STAR.TransposeFrom( L10 );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), X1_MC_STAR, L10Trans_MR_STAR, T(1), X0 );

        L11_STAR_STAR = L11;
        X1_VC_STAR = X1_MC_STAR;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
}

// Right Lower Normal (Non)Unit Trmm
//   X := X tril(L), and
//   X := X trilu(L)
template<typename T>
inline void
TrmmRLN
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("internal::TrmmRLN"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        TrmmRLNA( diag, alpha, L, X );
    else
        TrmmRLNC( diag, alpha, L, X );
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_BLAS_TRMM_RLN_HPP
