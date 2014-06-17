/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/



namespace El {
namespace symm {

template<typename T>
void LocalAccumulateRU
( Orientation orientation, T alpha,
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,STAR,MC  >& B_STAR_MC,
  const DistMatrix<T,MR,  STAR>& BTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZTrans_MC_STAR,
        DistMatrix<T,MR,  STAR>& ZTrans_MR_STAR )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::LocalAccumulateRU");
        if( A.Grid() != B_STAR_MC.Grid() ||
            B_STAR_MC.Grid() != BTrans_MR_STAR.Grid() ||
            BTrans_MR_STAR.Grid() != ZTrans_MC_STAR.Grid() ||
            ZTrans_MC_STAR.Grid() != ZTrans_MR_STAR.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
        if( A.Height() != A.Width() ||
            A.Height() != B_STAR_MC.Width() ||
            A.Height() != BTrans_MR_STAR.Height() ||
            A.Height() != ZTrans_MC_STAR.Height() ||
            A.Height() != ZTrans_MR_STAR.Height() ||
            B_STAR_MC.Height() != BTrans_MR_STAR.Width() ||
            BTrans_MR_STAR.Width() != ZTrans_MC_STAR.Width() ||
            ZTrans_MC_STAR.Width() != ZTrans_MR_STAR.Width() )
            LogicError
            ("Nonconformal:\n",
             DimsString(A,"A"),"\n",
             DimsString(B_STAR_MC,"B[* ,MC]"),"\n",
             DimsString(BTrans_MR_STAR,"B'[MR,* ]"),"\n",
             DimsString(ZTrans_MC_STAR,"Z'[MC,* ]"),"\n",
             DimsString(ZTrans_MR_STAR,"Z'[MR,* ]"));
        if( B_STAR_MC.RowAlign() != A.ColAlign() ||
            BTrans_MR_STAR.ColAlign() != A.RowAlign() ||
            ZTrans_MC_STAR.ColAlign() != A.ColAlign() ||
            ZTrans_MR_STAR.ColAlign() != A.RowAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MC>
        BL_STAR_MC(g), BR_STAR_MC(g),
        B0_STAR_MC(g), B1_STAR_MC(g), B2_STAR_MC(g);

    DistMatrix<T,MR,STAR>
        BTTrans_MR_STAR(g),  B0Trans_MR_STAR(g),
        BBTrans_MR_STAR(g),  B1Trans_MR_STAR(g),
                             B2Trans_MR_STAR(g);

    DistMatrix<T,MC,STAR>
        ZTTrans_MC_STAR(g),  Z0Trans_MC_STAR(g),
        ZBTrans_MC_STAR(g),  Z1Trans_MC_STAR(g),
                             Z2Trans_MC_STAR(g);

    DistMatrix<T,MR,STAR>
        ZBTrans_MR_STAR(g),  Z0Trans_MR_STAR(g),
        ZTTrans_MR_STAR(g),  Z1Trans_MR_STAR(g),
                             Z2Trans_MR_STAR(g);

    const Int ratio = Max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( B_STAR_MC,  BL_STAR_MC, BR_STAR_MC, 0 );
    LockedPartitionDown
    ( BTrans_MR_STAR, BTTrans_MR_STAR,
                      BBTrans_MR_STAR, 0 );
    PartitionDown
    ( ZTrans_MC_STAR, ZTTrans_MC_STAR,
                      ZBTrans_MC_STAR, 0 );
    PartitionDown
    ( ZTrans_MR_STAR, ZTTrans_MR_STAR,
                      ZBTrans_MR_STAR, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight
        ( BL_STAR_MC, /**/ BR_STAR_MC,
          B0_STAR_MC, /**/ B1_STAR_MC, B2_STAR_MC );

        LockedRepartitionDown
        ( BTTrans_MR_STAR,  B0Trans_MR_STAR,
         /***************/ /***************/
                            B1Trans_MR_STAR,
          BBTrans_MR_STAR,  B2Trans_MR_STAR );

        RepartitionDown
        ( ZTTrans_MC_STAR,  Z0Trans_MC_STAR,
         /***************/ /***************/
                            Z1Trans_MC_STAR,
          ZBTrans_MC_STAR,  Z2Trans_MC_STAR );

        RepartitionDown
        ( ZTTrans_MR_STAR,  Z0Trans_MR_STAR,
         /***************/ /***************/
                            Z1Trans_MR_STAR,
          ZBTrans_MR_STAR,  Z2Trans_MR_STAR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        MakeTriangular( UPPER, D11 );
        LocalGemm
        ( orientation, orientation,
          alpha, D11, B1_STAR_MC, T(1), Z1Trans_MR_STAR );
        SetDiagonal( D11, T(0) );

        LocalGemm
        ( NORMAL, NORMAL, alpha, D11, B1Trans_MR_STAR, T(1), Z1Trans_MC_STAR );

        LocalGemm
        ( orientation, orientation, 
          alpha, A12, B1_STAR_MC, T(1), Z2Trans_MR_STAR );

        LocalGemm
        ( NORMAL, NORMAL, alpha, A12, B2Trans_MR_STAR, T(1), Z1Trans_MC_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight
        ( BL_STAR_MC,             /**/ BR_STAR_MC,
          B0_STAR_MC, B1_STAR_MC, /**/ B2_STAR_MC );

        SlideLockedPartitionDown
        ( BTTrans_MR_STAR,  B0Trans_MR_STAR,
                            B1Trans_MR_STAR,
         /***************/ /***************/
          BBTrans_MR_STAR,  B2Trans_MR_STAR );

        SlidePartitionDown
        ( ZTTrans_MC_STAR,  Z0Trans_MC_STAR,
                            Z1Trans_MC_STAR,
         /***************/ /***************/
          ZBTrans_MC_STAR,  Z2Trans_MC_STAR );

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
RUA
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::RUA");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<T>
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);
    DistMatrix<T>
        CT(g),  C0(g),
        CB(g),  C1(g),
                C2(g);

    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,VC,  STAR> B1Trans_VC_STAR(g);
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);
    DistMatrix<T,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<T,MC,  MR  > Z1Trans(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    B1Trans_MR_STAR.AlignWith( A );
    B1Trans_VC_STAR.AlignWith( A );
    B1_STAR_MC.AlignWith( A );
    Z1Trans_MC_STAR.AlignWith( A );
    Z1Trans_MR_STAR.AlignWith( A );

    Matrix<T> Z1Local;

    Scale( beta, C );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( CT.Height() < C.Height() )
    {
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

        Z1Trans_MR_MC.AlignWith( C1 );
        //--------------------------------------------------------------------//
        B1.TransposeColAllGather( B1Trans_MR_STAR, conjugate );
        B1Trans_VC_STAR = B1Trans_MR_STAR;
        B1Trans_VC_STAR.TransposePartialColAllGather( B1_STAR_MC, conjugate );
        Zeros( Z1Trans_MC_STAR, C1.Width(), C1.Height() );
        Zeros( Z1Trans_MR_STAR, C1.Width(), C1.Height() );
        LocalAccumulateRU
        ( orientation, alpha, A, B1_STAR_MC, B1Trans_MR_STAR, 
          Z1Trans_MC_STAR, Z1Trans_MR_STAR );

        Z1Trans.RowSumScatterFrom( Z1Trans_MC_STAR );
        Z1Trans_MR_MC = Z1Trans;
        Z1Trans_MR_MC.RowSumScatterUpdate( T(1), Z1Trans_MR_STAR );
        Transpose( Z1Trans_MR_MC.LockedMatrix(), Z1Local, conjugate );
        Axpy( T(1), Z1Local, C1.Matrix() );
        //--------------------------------------------------------------------//

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
RUC
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("symm::RUC");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must be distributed on the same grid");
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Matrix views
    DistMatrix<T> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AB1(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  A1R(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);
    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g),
                  CLeft(g), CRight(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> AB1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > AB1Trans_STAR_MR(g);
    DistMatrix<T,MR,  STAR> A1RTrans_MR_STAR(g);

    B1_MC_STAR.AlignWith( C );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    PartitionRight( C, CL, CR, 0 );
    while( CR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );

        LockedView1x2( A1R, A11, A12 );
        LockedView2x1( AB1, A01, A11 );

        View1x2( CLeft, C0, C1 );
        View1x2( CRight, C1, C2 );

        AB1_VR_STAR.AlignWith( CLeft );
        AB1Trans_STAR_MR.AlignWith( CLeft );
        A1RTrans_MR_STAR.AlignWith( CRight );
        //--------------------------------------------------------------------//
        B1_MC_STAR = B1;

        AB1_VR_STAR = AB1;
        AB1_VR_STAR.TransposePartialColAllGather( AB1Trans_STAR_MR, conjugate );
        A1R.TransposeColAllGather( A1RTrans_MR_STAR, conjugate );
        MakeTriangular( LOWER, A1RTrans_MR_STAR );
        MakeTrapezoidal
        ( LOWER, AB1Trans_STAR_MR, 
          AB1Trans_STAR_MR.Width()-AB1Trans_STAR_MR.Height()-1 );

        LocalGemm
        ( NORMAL, orientation, 
          alpha, B1_MC_STAR, A1RTrans_MR_STAR, T(1), CRight );

        LocalGemm
        ( NORMAL, NORMAL,
          alpha, B1_MC_STAR, AB1Trans_STAR_MR, T(1), CLeft );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

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
RU
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("symm::RU"))
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Height() )
        symm::RUA( alpha, A, B, beta, C, conjugate );
    else
        symm::RUC( alpha, A, B, beta, C, conjugate );
}

} // namespace symm
} // namespace El
