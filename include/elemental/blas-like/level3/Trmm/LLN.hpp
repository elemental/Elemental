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
namespace internal {

template<typename T>
inline void
LocalTrmmAccumulateLLN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& L,
  const DistMatrix<T,STAR,MR  >& XTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::LocalTrmmAccumulateLLN");
        if( L.Grid() != XTrans_STAR_MR.Grid() ||
            XTrans_STAR_MR.Grid() != Z_MC_STAR.Grid() )
            LogicError("{L,X,Z} must be distributed over the same grid");
        if( L.Height() != L.Width() ||
            L.Height() != XTrans_STAR_MR.Width() ||
            L.Height() != Z_MC_STAR.Height() ||
            XTrans_STAR_MR.Height() != Z_MC_STAR.Width() )
            LogicError
            ("Nonconformal LocalTrmmAccumulateLLN: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X^H/T[* ,MR] ~ ",XTrans_STAR_MR.Height()," x ",
                                 XTrans_STAR_MR.Width(),"\n",
             "  Z[MC,* ] ~ ",Z_MC_STAR.Height()," x ",Z_MC_STAR.Width());
        if( XTrans_STAR_MR.RowAlign() != L.RowAlign() ||
            Z_MC_STAR.ColAlign() != L.ColAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

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
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionRight
    ( XTrans_STAR_MR, XLTrans_STAR_MR, XRTrans_STAR_MR, 0 );
    PartitionDown
    ( Z_MC_STAR, ZT_MC_STAR,
                 ZB_MC_STAR, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionRight
        ( XLTrans_STAR_MR, /**/ XRTrans_STAR_MR,
          X0Trans_STAR_MR, /**/ X1Trans_STAR_MR, 
                                X2Trans_STAR_MR );

        RepartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
         /**********/ /**********/
                       Z1_MC_STAR,
          ZB_MC_STAR,  Z2_MC_STAR );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        MakeTriangular( LOWER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm
        ( NORMAL, orientation, alpha, D11, X1Trans_STAR_MR, T(1), Z1_MC_STAR );
        LocalGemm
        ( NORMAL, orientation, alpha, L21, X1Trans_STAR_MR, T(1), Z2_MC_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

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
TrmmLLNA
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmLLNA");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() != L.Width() || L.Width() != X.Height() )
            LogicError
            ("Nonconformal TrmmLLNA: \n"
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Grid& g = L.Grid();

    DistMatrix<T>
        XL(g), XR(g),
        X0(g), X1(g), X2(g);

    DistMatrix<T,VR,  STAR> X1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > X1Trans_STAR_MR(g);
    DistMatrix<T,MC,  STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( L );
    X1Trans_STAR_MR.AlignWith( L );
    Z1_MC_STAR.AlignWith( L );

    PartitionRight( X, XL, XR, 0 );
    while( XL.Width() < X.Width() )
    {
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        //--------------------------------------------------------------------//
        X1_VR_STAR = X1;
        X1_VR_STAR.TransposePartialColAllGather( X1Trans_STAR_MR );
        Zeros( Z1_MC_STAR, X1.Height(), X1.Width() );
        LocalTrmmAccumulateLLN
        ( TRANSPOSE, diag, alpha, L, X1Trans_STAR_MR, Z1_MC_STAR );
        X1.RowSumScatterFrom( Z1_MC_STAR );
        //--------------------------------------------------------------------//

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
}

template<typename T>
inline void
TrmmLLNCOld
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmLLNCOld");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() != L.Width() || L.Width() != X.Height() )
            LogicError
            ("Nonconformal TrmmLLNC: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);
    DistMatrix<T> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        L10_STAR_MC.AlignWith( X0 );
        D1Trans_MR_STAR.AlignWith( X1 );
        D1Trans_MR_MC.AlignWith( X1 );
        D1.AlignWith( X1 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11;
        X1_STAR_VR = X1;
        LocalTrmm( LEFT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;

        L10_STAR_MC = L10;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, T(1), X0, L10_STAR_MC, D1Trans_MR_STAR );
        D1Trans_MR_MC.RowSumScatterFrom( D1Trans_MR_STAR );
        Zeros( D1, X1.Height(), X1.Width() );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix() );
        Axpy( T(1), D1, X1 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12, 
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
}

template<typename T>
inline void
TrmmLLNC
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrmmLLNC");
        if( L.Grid() != X.Grid() )
            LogicError("L and X must be distributed over the same grid");
        if( L.Height() != L.Width() || L.Width() != X.Height() )
            LogicError
            ("Nonconformal TrmmLLNC: \n",
             "  L ~ ",L.Height()," x ",L.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width());
    )
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);
    DistMatrix<T> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        L21_MC_STAR.AlignWith( X2 );
        X1Trans_MR_STAR.AlignWith( X2 );
        X1_STAR_VR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        L21_MC_STAR = L21;
        X1.TransposeColAllGather( X1Trans_MR_STAR );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), L21_MC_STAR, X1Trans_MR_STAR, T(1), X2 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.TransposePartialRowFilterFrom( X1Trans_MR_STAR );
        LocalTrmm( LEFT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12, 
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
}

// Left Lower Normal (Non)Unit Trmm 
//   X := tril(L)  X, or
//   X := trilu(L) X
template<typename T>
inline void
TrmmLLN
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("internal::TrmmLLN"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        TrmmLLNA( diag, alpha, L, X );
    else
        TrmmLLNC( diag, alpha, L, X );
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRMM_LLN_HPP
