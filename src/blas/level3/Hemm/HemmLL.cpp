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
elemental::blas::internal::HemmLL
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmLL");
#endif
    blas::internal::HemmLLC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::HemmLLA
( T alpha, const DistMatrix<T,MC,MR>& A, 
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmLLA");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.GetGrid();

    DistMatrix<T,MC,MR> 
        BL(g), BR(g),
        B0(g), B1(g), B2(g);

    DistMatrix<T,MC,MR>
        CL(g), CR(g),
        C0(g), C1(g), C2(g);

    DistMatrix<T,MC,Star> B1_MC_Star(g);
    DistMatrix<T,MR,Star> B1_MR_Star(g);
    DistMatrix<T,MC,MR  > Z1(g);
    DistMatrix<T,MC,Star> Z1_MC_Star(g);
    DistMatrix<T,MR,Star> Z1_MR_Star(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    blas::Scal( beta, C );
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

        B1_MC_Star.AlignWith( A );
        B1_MR_Star.AlignWith( A );
        Z1_MC_Star.AlignWith( A );
        Z1_MR_Star.AlignWith( A );
        Z1.AlignWith( C1 );
        Z1_MC_Star.ResizeTo( C1.Height(), C1.Width() );
        Z1_MR_Star.ResizeTo( C1.Height(), C1.Width() );
        Z1_MC_Star.SetToZero();
        Z1_MR_Star.SetToZero();
        //--------------------------------------------------------------------//
        B1_MC_Star = B1;
        B1_MR_Star = B1_MC_Star;
        blas::internal::LocalHemmAccumulateLL
        ( alpha, A, B1_MC_Star, B1_MR_Star, Z1_MC_Star, Z1_MR_Star );

        Z1_MR_MC.SumScatterFrom( Z1_MR_Star );
        Z1 = Z1_MR_MC;
        Z1.SumScatterUpdate( (T)1, Z1_MC_Star );
        blas::Axpy( (T)1, Z1, C1 );
        //--------------------------------------------------------------------//
        B1_MC_Star.FreeAlignments();
        B1_MR_Star.FreeAlignments();
        Z1_MC_Star.FreeAlignments();
        Z1_MR_Star.FreeAlignments();
        Z1.FreeAlignments();

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
elemental::blas::internal::HemmLLC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmLLC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AColPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  ARowPan(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> 
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);

    DistMatrix<T,MC,MR> 
        CT(g),  C0(g),  CAbove(g),
        CB(g),  C1(g),  CBelow(g),
                C2(g);

    // Temporary distributions
    DistMatrix<T,MC,Star> AColPan_MC_Star(g);
    DistMatrix<T,Star,MC> ARowPan_Star_MC(g);
    DistMatrix<T,Star,MR> B1_Star_MR(g);

    // Start the algorithm
    blas::Scal( beta, C );
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

        ARowPan.LockedView1x2( A10, A11 );

        AColPan.LockedView2x1
        ( A11,
          A21 );

        CAbove.View2x1
        ( C0,
          C1 );

        CBelow.View2x1
        ( C1,
          C2 );

        AColPan_MC_Star.AlignWith( CBelow );
        ARowPan_Star_MC.AlignWith( CAbove );
        B1_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        AColPan_MC_Star = AColPan;
        ARowPan_Star_MC = ARowPan;
        AColPan_MC_Star.MakeTrapezoidal( Left, Lower );
        ARowPan_Star_MC.MakeTrapezoidal( Right, Lower, -1 );

        B1_Star_MR = B1;

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, AColPan_MC_Star, B1_Star_MR, (T)1, CBelow );

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          alpha, ARowPan_Star_MC, B1_Star_MR, (T)1, CAbove );
        //--------------------------------------------------------------------//
        AColPan_MC_Star.FreeAlignments();
        ARowPan_Star_MC.FreeAlignments();
        B1_Star_MR.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalHemmAccumulateLL
( T alpha,  
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,Star>& B_MC_Star,
  const DistMatrix<T,MR,Star>& B_MR_Star,
        DistMatrix<T,MC,Star>& Z_MC_Star,
        DistMatrix<T,MR,Star>& Z_MR_Star )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalHemmAccumulateLL");
    if( A.GetGrid() != B_MC_Star.GetGrid() || 
        B_MC_Star.GetGrid() != B_MR_Star.GetGrid() ||
        B_MR_Star.GetGrid() != Z_MC_Star.GetGrid() ||
        Z_MC_Star.GetGrid() != Z_MR_Star.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != A.Width() || 
        A.Height() != B_MC_Star.Height() ||
        A.Height() != B_MR_Star.Height() ||
        A.Height() != Z_MC_Star.Height() ||
        A.Height() != Z_MR_Star.Height() ||
        B_MC_Star.Width() != B_MR_Star.Width() ||
        B_MR_Star.Width() != Z_MC_Star.Width() ||
        Z_MC_Star.Width() != Z_MR_Star.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalHemmAccumulateLL: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B[MC,* ] ~ " << B_MC_Star.Height() << " x " 
                               << B_MC_Star.Width() << endl
            << "  B[MR,* ] ~ " << B_MR_Star.Height() << " x " 
                               << B_MR_Star.Width() << endl
            << "  Z[MC,* ] ~ " << Z_MC_Star.Height() << " x " 
                               << Z_MC_Star.Width() << endl
            << "  Z[MR,* ] ` " << Z_MR_Star.Height() << " x " 
                               << Z_MR_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( B_MC_Star.ColAlignment() != A.ColAlignment() || 
        B_MR_Star.ColAlignment() != A.RowAlignment() ||
        Z_MC_Star.ColAlignment() != A.ColAlignment() ||
        Z_MR_Star.ColAlignment() != A.RowAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,MC,Star> 
        BT_MC_Star(g),  B0_MC_Star(g),
        BB_MC_Star(g),  B1_MC_Star(g),
                        B2_MC_Star(g);

    DistMatrix<T,MR,Star>
        BT_MR_Star(g),  B0_MR_Star(g),
        BB_MR_Star(g),  B1_MR_Star(g),
                        B2_MR_Star(g);

    DistMatrix<T,MC,Star>
        ZT_MC_Star(g),  Z0_MC_Star(g),
        ZB_MC_Star(g),  Z1_MC_Star(g),
                        Z2_MC_Star(g);

    DistMatrix<T,MR,Star>
        ZT_MR_Star(g),  Z0_MR_Star(g),
        ZB_MR_Star(g),  Z1_MR_Star(g),
                        Z2_MR_Star(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDown
    ( B_MC_Star, BT_MC_Star,
                 BB_MC_Star, 0 );
    LockedPartitionDown
    ( B_MR_Star, BT_MR_Star,
                 BB_MR_Star, 0 );
    PartitionDown
    ( Z_MC_Star, ZT_MC_Star,
                 ZB_MC_Star, 0 );
    PartitionDown
    ( Z_MR_Star, ZT_MR_Star,
                 ZB_MR_Star, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown
        ( BT_MC_Star,  B0_MC_Star,
         /**********/ /**********/
                       B1_MC_Star,
          BB_MC_Star,  B2_MC_Star );

        LockedRepartitionDown
        ( BT_MR_Star,  B0_MR_Star,
         /**********/ /**********/
                       B1_MR_Star,
          BB_MR_Star,  B2_MR_Star );

        RepartitionDown
        ( ZT_MC_Star,  Z0_MC_Star,
         /**********/ /**********/
                       Z1_MC_Star,
          ZB_MC_Star,  Z2_MC_Star );

        RepartitionDown
        ( ZT_MR_Star,  Z0_MR_Star,
         /**********/ /**********/
                       Z1_MR_Star,
          ZB_MR_Star,  Z2_MR_Star );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Lower );
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, D11, B1_MR_Star, (T)1, Z1_MC_Star );
        D11.MakeTrapezoidal( Left, Lower, -1 );

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          alpha, D11, B1_MC_Star, (T)1, Z1_MR_Star );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A21, B1_MR_Star, (T)1, Z2_MC_Star );

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          alpha, A21, B2_MC_Star, (T)1, Z1_MR_Star );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown
        ( BT_MC_Star,  B0_MC_Star,
                       B1_MC_Star,
         /**********/ /**********/
          BB_MC_Star,  B2_MC_Star );

        SlideLockedPartitionDown
        ( BT_MR_Star,  B0_MR_Star,
                       B1_MR_Star,
         /**********/ /**********/
          BB_MR_Star,  B2_MR_Star );

        SlidePartitionDown
        ( ZT_MC_Star,  Z0_MC_Star,
                       Z1_MC_Star,
         /**********/ /**********/
          ZB_MC_Star,  Z2_MC_Star );

        SlidePartitionDown
        ( ZT_MR_Star,  Z0_MR_Star,
                       Z1_MR_Star,
         /**********/ /**********/
          ZB_MR_Star,  Z2_MR_Star );
    }

    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::HemmLL
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::HemmLLA
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::HemmLL
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::HemmLLA
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::HemmLL
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::HemmLLA
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::HemmLL
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::HemmLLA
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

