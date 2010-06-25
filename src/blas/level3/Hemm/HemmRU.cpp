/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::HemmRU
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmRU");
#endif
    blas::internal::HemmRUC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::HemmRUA
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmRUA");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.GetGrid();

    DistMatrix<T,MC,MR>
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);

    DistMatrix<T,MC,MR>
        CT(g),  C0(g),
        CB(g),  C1(g),
                C2(g);

    DistMatrix<T,MR,  Star> B1Herm_MR_Star(g);
    DistMatrix<T,VC,  Star> B1Herm_VC_Star(g);
    DistMatrix<T,Star,MC  > B1_Star_MC(g);
    DistMatrix<T,MC  ,MR  > Z1(g);
    DistMatrix<T,Star,MC  > Z1_Star_MC(g);
    DistMatrix<T,Star,MR  > Z1_Star_MR(g);
    DistMatrix<T,MR,  MC  > Z1_MR_MC(g);

    blas::Scal( beta, C );
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

        B1Herm_MR_Star.AlignWith( A );
        B1Herm_VC_Star.AlignWith( A );
        B1_Star_MC.AlignWith( A );
        Z1_Star_MC.AlignWith( A );
        Z1_Star_MR.AlignWith( A );
        Z1.AlignWith( C1 );
        Z1_Star_MC.ResizeTo( C1.Height(), C1.Width() );
        Z1_Star_MR.ResizeTo( C1.Height(), C1.Width() );
        Z1_Star_MC.SetToZero();
        Z1_Star_MR.SetToZero();
        //--------------------------------------------------------------------//
        B1Herm_MR_Star.ConjugateTransposeFrom( B1 );
        B1Herm_VC_Star = B1Herm_MR_Star;
        B1_Star_MC.ConjugateTransposeFrom( B1Herm_VC_Star );
        blas::internal::LocalHemmAccumulateRU
        ( alpha, A, B1_Star_MC, B1Herm_MR_Star, Z1_Star_MC, Z1_Star_MR );

        Z1_MR_MC.SumScatterFrom( Z1_Star_MC );
        Z1 = Z1_MR_MC;
        Z1.SumScatterUpdate( (T)1, Z1_Star_MR );
        blas::Axpy( (T)1, Z1, C1 );
        //--------------------------------------------------------------------//
        B1Herm_MR_Star.FreeAlignments();
        B1Herm_VC_Star.FreeAlignments();
        B1_Star_MC.FreeAlignments();
        Z1_Star_MC.FreeAlignments();
        Z1_Star_MR.FreeAlignments();
        Z1.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::HemmRUC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmRUC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed on the same grid." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AColPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  ARowPan(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    DistMatrix<T,MC,MR> CL(g), CR(g),
                        C0(g), C1(g), C2(g),
                        CLeft(g), CRight(g);

    // Temporary distributions
    DistMatrix<T,MC,Star> B1_MC_Star(g);
    DistMatrix<T,VR,  Star> AColPan_VR_Star(g);
    DistMatrix<T,Star,MR  > AColPanHerm_Star_MR(g);
    DistMatrix<T,MR,  Star> ARowPanHerm_MR_Star(g);

    // Start the algorithm
    blas::Scal( beta, C );
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

        ARowPan.LockedView1x2( A11, A12 );

        AColPan.LockedView2x1
        ( A01,
          A11 );

        CLeft.View1x2( C0, C1 );

        CRight.View1x2( C1, C2 );

        B1_MC_Star.AlignWith( C );
        AColPan_VR_Star.AlignWith( CLeft );
        AColPanHerm_Star_MR.AlignWith( CLeft );
        ARowPanHerm_MR_Star.AlignWith( CRight );
        //--------------------------------------------------------------------//
        B1_MC_Star = B1;

        AColPan_VR_Star = AColPan;
        AColPanHerm_Star_MR.ConjugateTransposeFrom( AColPan_VR_Star );
        ARowPanHerm_MR_Star.ConjugateTransposeFrom( ARowPan );
        ARowPanHerm_MR_Star.MakeTrapezoidal( Left, Lower );
        AColPanHerm_Star_MR.MakeTrapezoidal( Right, Lower, -1 );

        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          alpha, B1_MC_Star, ARowPanHerm_MR_Star, (T)1, CRight );

        blas::internal::LocalGemm
        ( Normal, Normal,
          alpha, B1_MC_Star, AColPanHerm_Star_MR, (T)1, CLeft );
        //--------------------------------------------------------------------//
        B1_MC_Star.FreeAlignments();
        AColPan_VR_Star.FreeAlignments();
        AColPanHerm_Star_MR.FreeAlignments();
        ARowPanHerm_MR_Star.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalHemmAccumulateRU
( T alpha,
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,Star,MC  >& B_Star_MC,
  const DistMatrix<T,MR,  Star>& BHerm_MR_Star,
        DistMatrix<T,Star,MC  >& Z_Star_MC,
        DistMatrix<T,Star,MR  >& Z_Star_MR )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalHemmAccumulateRU");
    if( A.GetGrid() != B_Star_MC.GetGrid() ||
        B_Star_MC.GetGrid() != BHerm_MR_Star.GetGrid() ||
        BHerm_MR_Star.GetGrid() != Z_Star_MC.GetGrid() ||
        Z_Star_MC.GetGrid() != Z_Star_MR.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != A.Width() ||
        A.Height() != B_Star_MC.Width() ||
        A.Height() != BHerm_MR_Star.Height() ||
        A.Height() != Z_Star_MC.Width() ||
        A.Height() != Z_Star_MR.Width() ||
        B_Star_MC.Height() != BHerm_MR_Star.Width() ||
        BHerm_MR_Star.Width() != Z_Star_MC.Height() ||
        Z_Star_MC.Height() != Z_Star_MR.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalHemmAccumulateRU: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B[* ,MC] ~ " << B_Star_MC.Height() << " x "
                               << B_Star_MC.Width() << endl
            << "  B^H[MR,* ] ~ " << BHerm_MR_Star.Height() << " x "
                                 << BHerm_MR_Star.Width() << endl
            << "  Z[* ,MC] ~ " << Z_Star_MC.Height() << " x "
                               << Z_Star_MC.Width() << endl
            << "  Z[* ,MR] ~ " << Z_Star_MR.Height() << " x "
                               << Z_Star_MR.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( B_Star_MC.RowAlignment() != A.ColAlignment() ||
        BHerm_MR_Star.ColAlignment() != A.RowAlignment() ||
        Z_Star_MC.RowAlignment() != A.ColAlignment() ||
        Z_Star_MR.RowAlignment() != A.RowAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,Star,MC>
        BL_Star_MC(g), BR_Star_MC(g),
        B0_Star_MC(g), B1_Star_MC(g), B2_Star_MC(g);

    DistMatrix<T,MR,Star>
        BHermT_MR_Star(g),  BHerm0_MR_Star(g),
        BHermB_MR_Star(g),  BHerm1_MR_Star(g),
                            BHerm2_MR_Star(g);

    DistMatrix<T,Star,MC>
        ZL_Star_MC(g), ZR_Star_MC(g),
        Z0_Star_MC(g), Z1_Star_MC(g), Z2_Star_MC(g);

    DistMatrix<T,Star,MR>
        ZL_Star_MR(g), ZR_Star_MR(g),
        Z0_Star_MR(g), Z1_Star_MR(g), Z2_Star_MR(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( B_Star_MC,  BL_Star_MC, BR_Star_MC, 0 );
    LockedPartitionDown
    ( BHerm_MR_Star, BHermT_MR_Star,
                     BHermB_MR_Star, 0 );
    PartitionRight( Z_Star_MC,  ZL_Star_MC, ZR_Star_MC, 0 );
    PartitionRight( Z_Star_MR,  ZL_Star_MR, ZR_Star_MR, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight
        ( BL_Star_MC, /**/ BR_Star_MC,
          B0_Star_MC, /**/ B1_Star_MC, B2_Star_MC );

        LockedRepartitionDown
        ( BHermT_MR_Star,  BHerm0_MR_Star,
         /**************/ /**************/
                           BHerm1_MR_Star,
          BHermB_MR_Star,  BHerm2_MR_Star );

        RepartitionRight
        ( ZL_Star_MC, /**/ ZR_Star_MC,
          Z0_Star_MC, /**/ Z1_Star_MC, Z2_Star_MC );

        RepartitionRight
        ( ZL_Star_MR, /**/ ZR_Star_MR,
          Z0_Star_MR, /**/ Z1_Star_MR, Z2_Star_MR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Upper );
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, B1_Star_MC, D11, (T)1, Z1_Star_MR );
        D11.MakeTrapezoidal( Left, Upper, 1 );

        blas::internal::LocalGemm
        ( ConjugateTranspose, ConjugateTranspose,
          alpha, BHerm1_MR_Star, D11, (T)1, Z1_Star_MC );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, B1_Star_MC, A12, (T)1, Z2_Star_MR );

        blas::internal::LocalGemm
        ( ConjugateTranspose, ConjugateTranspose,
          alpha, BHerm2_MR_Star, A12, (T)1, Z1_Star_MC );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight
        ( BL_Star_MC,             /**/ BR_Star_MC,
          B0_Star_MC, B1_Star_MC, /**/ B2_Star_MC );

        SlideLockedPartitionDown
        ( BHermT_MR_Star,  BHerm0_MR_Star,
                           BHerm1_MR_Star,
         /**************/ /**************/
          BHermB_MR_Star,  BHerm2_MR_Star );

        SlidePartitionRight
        ( ZL_Star_MC,             /**/ ZR_Star_MC,
          Z0_Star_MC, Z1_Star_MC, /**/ Z2_Star_MC );

        SlidePartitionRight
        ( ZL_Star_MR,             /**/ ZR_Star_MR,
          Z0_Star_MR, Z1_Star_MR, /**/ Z2_Star_MR );
    }

    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::HemmRU
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::HemmRUA
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::HemmRU
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::HemmRUA
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::HemmRU
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::HemmRUA
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::HemmRU
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::HemmRUA
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

