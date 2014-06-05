/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include EL_ZEROS_INC

namespace El {
namespace symm {

template<typename T>
void LocalAccumulateLU
( Orientation orientation, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,  STAR>& B_MC_STAR,
  const DistMatrix<T,STAR,MR  >& BTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::LocalAccumulateLU");
        if( A.Grid() != B_MC_STAR.Grid() ||
            B_MC_STAR.Grid() != BTrans_STAR_MR.Grid() ||
            BTrans_STAR_MR.Grid() != Z_MC_STAR.Grid() ||
            Z_MC_STAR.Grid() != Z_MR_STAR.Grid() )
            LogicError("{A,B,Z} must be distributed over the same grid");
        if( A.Height() != A.Width() ||
            A.Height() != B_MC_STAR.Height() ||
            A.Height() != BTrans_STAR_MR.Width() ||
            A.Height() != Z_MC_STAR.Height() ||
            A.Height() != Z_MR_STAR.Height() ||
            B_MC_STAR.Width() != BTrans_STAR_MR.Height() ||
            BTrans_STAR_MR.Height() != Z_MC_STAR.Width() ||
            Z_MC_STAR.Width() != Z_MR_STAR.Width() )
            LogicError
            ("Nonconformal:\n",
             DimsString(A,"A"),"\n",
             DimsString(B_MC_STAR,"B[MC,* ]"),"\n",
             DimsString(BTrans_STAR_MR,"B'[* ,MR]"),"\n",
             DimsString(Z_MC_STAR,"Z[MC,* ]"),"\n",
             DimsString(Z_MR_STAR,"Z[MR,* ]"));
        if( B_MC_STAR.ColAlign() != A.ColAlign() ||
            BTrans_STAR_MR.RowAlign() != A.RowAlign() ||
            Z_MC_STAR.ColAlign() != A.ColAlign() ||
            Z_MR_STAR.ColAlign() != A.RowAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = A.Grid();

    DistMatrix<T>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,MC,STAR>
        BT_MC_STAR(g),  B0_MC_STAR(g),
        BB_MC_STAR(g),  B1_MC_STAR(g),
                        B2_MC_STAR(g);

    DistMatrix<T,STAR,MR>
        BLTrans_STAR_MR(g), BRTrans_STAR_MR(g),
        B0Trans_STAR_MR(g), B1Trans_STAR_MR(g),
        B2Trans_STAR_MR(g);

    DistMatrix<T,MC,STAR>
        ZT_MC_STAR(g),  Z0_MC_STAR(g),
        ZB_MC_STAR(g),  Z1_MC_STAR(g),
                        Z2_MC_STAR(g);

    DistMatrix<T,MR,STAR>
        ZT_MR_STAR(g),  Z0_MR_STAR(g),
        ZB_MR_STAR(g),  Z1_MR_STAR(g),
                        Z2_MR_STAR(g);

    const Int ratio = Max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDown
    ( B_MC_STAR, BT_MC_STAR,
                 BB_MC_STAR, 0 );
    LockedPartitionRight
    ( BTrans_STAR_MR, BLTrans_STAR_MR, BRTrans_STAR_MR, 0 );
    PartitionDown
    ( Z_MC_STAR, ZT_MC_STAR,
                 ZB_MC_STAR, 0 );
    PartitionDown
    ( Z_MR_STAR, ZT_MR_STAR,
                 ZB_MR_STAR, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
          /************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown
        ( BT_MC_STAR,  B0_MC_STAR,
         /**********/ /**********/
                       B1_MC_STAR,
          BB_MC_STAR,  B2_MC_STAR );

        LockedRepartitionRight
        ( BLTrans_STAR_MR, /**/ BRTrans_STAR_MR,
          B0Trans_STAR_MR, /**/ B1Trans_STAR_MR,
                                B2Trans_STAR_MR );

        RepartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
         /**********/ /**********/
                       Z1_MC_STAR,
          ZB_MC_STAR,  Z2_MC_STAR );

        RepartitionDown
        ( ZT_MR_STAR,  Z0_MR_STAR,
         /**********/ /**********/
                       Z1_MR_STAR,
          ZB_MR_STAR,  Z2_MR_STAR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        MakeTriangular( UPPER, D11 );
        LocalGemm
        ( NORMAL, orientation, alpha, D11, B1Trans_STAR_MR, T(1), Z1_MC_STAR );
        SetDiagonal( D11, T(0) );

        LocalGemm
        ( orientation, NORMAL, alpha, D11, B1_MC_STAR, T(1), Z1_MR_STAR );

        LocalGemm
        ( NORMAL, orientation, alpha, A12, B2Trans_STAR_MR, T(1), Z1_MC_STAR );

        LocalGemm
        ( orientation, NORMAL, alpha, A12, B1_MC_STAR, T(1), Z2_MR_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown
        ( BT_MC_STAR,  B0_MC_STAR,
                       B1_MC_STAR,
         /**********/ /**********/
          BB_MC_STAR,  B2_MC_STAR );

        SlideLockedPartitionRight
        ( BLTrans_STAR_MR,                  /**/ BRTrans_STAR_MR,
          B0Trans_STAR_MR, B1Trans_STAR_MR, /**/ B2Trans_STAR_MR );

        SlidePartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
                       Z1_MC_STAR,
         /**********/ /**********/
          ZB_MC_STAR,  Z2_MC_STAR );

        SlidePartitionDown
        ( ZT_MR_STAR,  Z0_MR_STAR,
                       Z1_MR_STAR,
         /**********/ /**********/
          ZB_MR_STAR,  Z2_MR_STAR );
    }
    PopBlocksizeStack();
}

template<typename T>
inline void
LUA
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::LUA");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<T>
        BL(g), BR(g),
        B0(g), B1(g), B2(g);

    DistMatrix<T>
        CL(g), CR(g),
        C0(g), C1(g), C2(g);

    DistMatrix<T,MC,STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);
    DistMatrix<T> Z1(g);
    DistMatrix<T,MC,STAR> Z1_MC_STAR(g);
    DistMatrix<T,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    B1_MC_STAR.AlignWith( A );
    B1_VR_STAR.AlignWith( A );
    B1Trans_STAR_MR.AlignWith( A );
    Z1_MC_STAR.AlignWith( A );
    Z1_MR_STAR.AlignWith( A );

    Scale( beta, C );
    LockedPartitionRight
    ( B, BL, BR, 0 );
    PartitionRight
    ( C, CL, CR, 0 );
    while( CL.Width() < C.Width() )
    {
        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );

        Z1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        B1_MC_STAR = B1;
        B1_VR_STAR = B1_MC_STAR;
        B1_VR_STAR.TransposePartialColAllGather( B1Trans_STAR_MR, conjugate );
        Zeros( Z1_MC_STAR, C1.Height(), C1.Width() );
        Zeros( Z1_MR_STAR, C1.Height(), C1.Width() );
        LocalAccumulateLU
        ( orientation,
          alpha, A, B1_MC_STAR, B1Trans_STAR_MR, Z1_MC_STAR, Z1_MR_STAR );

        Z1_MR_MC.RowSumScatterFrom( Z1_MR_STAR );
        Z1 = Z1_MR_MC;
        Z1.RowSumScatterUpdate( T(1), Z1_MC_STAR );
        Axpy( T(1), Z1, C1 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
}

template<typename T>
inline void
LUC
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::LUC");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Matrix views
    DistMatrix<T> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AB1(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  A1R(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);
    DistMatrix<T> CT(g),  C0(g),  CAbove(g),
                  CB(g),  C1(g),  CBelow(g),
                          C2(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> AB1_MC_STAR(g);
    DistMatrix<T,STAR,MC  > A1R_STAR_MC(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);

    B1Trans_MR_STAR.AlignWith( C );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( CB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );

        RepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );

        LockedView1x2( A1R, A11, A12 );
        LockedView2x1( AB1, A01, A11 );

        View2x1( CAbove, C0, C1 );
        View2x1( CBelow, C1, C2 );

        AB1_MC_STAR.AlignWith( CAbove );
        A1R_STAR_MC.AlignWith( CBelow );
        //--------------------------------------------------------------------//
        AB1_MC_STAR = AB1;
        A1R_STAR_MC = A1R;
        MakeTrapezoidal
        ( UPPER, AB1_MC_STAR, AB1_MC_STAR.Width()-AB1_MC_STAR.Height() );
        MakeTrapezoidal( UPPER, A1R_STAR_MC, 1 );

        B1.TransposeColAllGather( B1Trans_MR_STAR );

        LocalGemm
        ( NORMAL, TRANSPOSE, 
          alpha, AB1_MC_STAR, B1Trans_MR_STAR, T(1), CAbove );

        LocalGemm
        ( orientation, TRANSPOSE, 
          alpha, A1R_STAR_MC, B1Trans_MR_STAR, T(1), CBelow );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );

        SlidePartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
    }
}

template<typename T>
inline void
LU
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("symm::LU"))
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Width() )
        symm::LUA( alpha, A, B, beta, C, conjugate );
    else
        symm::LUC( alpha, A, B, beta, C, conjugate );
}

} // namespace symm
} // namespace El
