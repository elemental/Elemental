/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::HemmRL
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmRL");
#endif
    blas::internal::HemmRLC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::HemmRLA
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmRLA");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.Grid();

    DistMatrix<T,MC,MR>
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);

    DistMatrix<T,MC,MR>
        CT(g),  C0(g),
        CB(g),  C1(g),
                C2(g);

    DistMatrix<T,MR,  Star> BHerm1_MR_Star(g);
    DistMatrix<T,VC,  Star> BHerm1_VC_Star(g);
    DistMatrix<T,Star,MC  > B1_Star_MC(g);
    DistMatrix<T,MC,  Star> Z1Herm_MC_Star(g);
    DistMatrix<T,MR,  Star> Z1Herm_MR_Star(g);
    DistMatrix<T,MC,  MR  > Z1Herm(g);
    DistMatrix<T,MR,  MC  > Z1Herm_MR_MC(g);

    Matrix<T> Z1Local;

    blas::Scal( beta, C );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( BT.Height() < B.Height() )
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

        BHerm1_MR_Star.AlignWith( A );
        BHerm1_VC_Star.AlignWith( A );
        B1_Star_MC.AlignWith( A );
        Z1Herm_MC_Star.AlignWith( A );
        Z1Herm_MR_Star.AlignWith( A );
        Z1Herm_MR_MC.AlignWith( C1 );
        Z1Herm_MC_Star.ResizeTo( C1.Width(), C1.Height() );
        Z1Herm_MR_Star.ResizeTo( C1.Width(), C1.Height() );
        Z1Herm_MC_Star.SetToZero();
        Z1Herm_MR_Star.SetToZero();
        //--------------------------------------------------------------------//
        BHerm1_MR_Star.ConjugateTransposeFrom( B1 );
        BHerm1_VC_Star = BHerm1_MR_Star;
        B1_Star_MC.ConjugateTransposeFrom( BHerm1_VC_Star );
        blas::internal::LocalHemmAccumulateRL
        ( alpha, A, B1_Star_MC, BHerm1_MR_Star, 
          Z1Herm_MC_Star, Z1Herm_MR_Star );

        Z1Herm.SumScatterFrom( Z1Herm_MC_Star );
        Z1Herm_MR_MC = Z1Herm;
        Z1Herm_MR_MC.SumScatterUpdate( (T)1, Z1Herm_MR_Star );
        blas::ConjTrans( Z1Herm_MR_MC.LockedLocalMatrix(), Z1Local );
        blas::Axpy( (T)1, Z1Local, C1.LocalMatrix() );
        //--------------------------------------------------------------------//
        BHerm1_MR_Star.FreeAlignments();
        BHerm1_VC_Star.FreeAlignments();
        B1_Star_MC.FreeAlignments();
        Z1Herm_MC_Star.FreeAlignments();
        Z1Herm_MR_Star.FreeAlignments();
        Z1Herm_MR_MC.FreeAlignments();

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
elemental::blas::internal::HemmRLC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmRLC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.Grid();

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
    DistMatrix<T,MC,  Star> B1_MC_Star(g);
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

        ARowPan.LockedView1x2( A10, A11 );

        AColPan.LockedView2x1
        ( A11,
          A21 );

        CLeft.View1x2( C0, C1 );

        CRight.View1x2( C1, C2 );

        B1_MC_Star.AlignWith( C );
        AColPan_VR_Star.AlignWith( CRight );
        AColPanHerm_Star_MR.AlignWith( CRight );
        ARowPanHerm_MR_Star.AlignWith( CLeft );
        //--------------------------------------------------------------------//
        B1_MC_Star = B1;

        ARowPanHerm_MR_Star.ConjugateTransposeFrom( ARowPan );
        AColPan_VR_Star = AColPan;
        AColPanHerm_Star_MR.ConjugateTransposeFrom( AColPan_VR_Star );
        ARowPanHerm_MR_Star.MakeTrapezoidal( Right, Upper );
        AColPanHerm_Star_MR.MakeTrapezoidal( Left, Upper, 1 );

        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          alpha, B1_MC_Star, ARowPanHerm_MR_Star, (T)1, CLeft );

        blas::internal::LocalGemm
        ( Normal, Normal, 
          alpha, B1_MC_Star, AColPanHerm_Star_MR, (T)1, CRight );
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
elemental::blas::internal::LocalHemmAccumulateRL
( T alpha,
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,Star,MC  >& B_Star_MC,
  const DistMatrix<T,MR,  Star>& BHerm_MR_Star,
        DistMatrix<T,MC,  Star>& ZHerm_MC_Star,
        DistMatrix<T,MR,  Star>& ZHerm_MR_Star )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalHemmAccumulateRL");
    if( A.Grid() != B_Star_MC.Grid() ||
        B_Star_MC.Grid() != BHerm_MR_Star.Grid() ||
        BHerm_MR_Star.Grid() != ZHerm_MC_Star.Grid() ||
        ZHerm_MC_Star.Grid() != ZHerm_MR_Star.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != A.Width() ||
        A.Height() != B_Star_MC.Width() ||
        A.Height() != BHerm_MR_Star.Height() ||
        A.Height() != ZHerm_MC_Star.Height() ||
        A.Height() != ZHerm_MR_Star.Height() ||
        B_Star_MC.Height() != BHerm_MR_Star.Width() ||
        BHerm_MR_Star.Width() != ZHerm_MC_Star.Width() ||
        ZHerm_MC_Star.Width() != ZHerm_MR_Star.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalHemmAccumulateRL: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B[* ,MC] ~ " << B_Star_MC.Height() << " x " 
                               << B_Star_MC.Width() << endl
            << "  B^H[MR,* ] ~ " << BHerm_MR_Star.Height() << " x " 
                               << BHerm_MR_Star.Width() << endl
            << "  Z^H[MC,* ] ~ " << ZHerm_MC_Star.Height() << " x " 
                                 << ZHerm_MC_Star.Width() << endl
            << "  Z^H[MR,* ] ~ " << ZHerm_MR_Star.Height() << " x " 
                                 << ZHerm_MR_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( B_Star_MC.RowAlignment() != A.ColAlignment() ||
        BHerm_MR_Star.ColAlignment() != A.RowAlignment() ||
        ZHerm_MC_Star.ColAlignment() != A.ColAlignment() ||
        ZHerm_MR_Star.ColAlignment() != A.RowAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = A.Grid();

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

    DistMatrix<T,MC,Star>
        ZHermT_MC_Star(g),  ZHerm0_MC_Star(g),
        ZHermB_MC_Star(g),  ZHerm1_MC_Star(g),
                            ZHerm2_MC_Star(g);

    DistMatrix<T,MR,Star>
        ZHermT_MR_Star(g),  ZHerm0_MR_Star(g),
        ZHermB_MR_Star(g),  ZHerm1_MR_Star(g),
                            ZHerm2_MR_Star(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( B_Star_MC,  BL_Star_MC, BR_Star_MC, 0 );
    LockedPartitionDown
    ( BHerm_MR_Star, BHermT_MR_Star,
                     BHermB_MR_Star, 0 );
    PartitionDown
    ( ZHerm_MC_Star, ZHermT_MC_Star,
                     ZHermB_MC_Star, 0 );
    PartitionDown
    ( ZHerm_MR_Star, ZHermT_MR_Star,
                     ZHermB_MR_Star, 0 );
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

        RepartitionDown
        ( ZHermT_MC_Star,  ZHerm0_MC_Star,
         /**************/ /**************/
                           ZHerm1_MC_Star,
          ZHermB_MC_Star,  ZHerm2_MC_Star );

        RepartitionDown
        ( ZHermT_MR_Star,  ZHerm0_MR_Star,
         /**************/ /**************/
                           ZHerm1_MR_Star,
          ZHermB_MR_Star,  ZHerm2_MR_Star );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Lower );
        blas::internal::LocalGemm
        ( ConjugateTranspose, ConjugateTranspose, 
          alpha, D11, B1_Star_MC, (T)1, ZHerm1_MR_Star );
        D11.MakeTrapezoidal( Left, Lower, -1 );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, D11, BHerm1_MR_Star, (T)1, ZHerm1_MC_Star );

        blas::internal::LocalGemm
        ( ConjugateTranspose, ConjugateTranspose,
          alpha, A21, B2_Star_MC, (T)1, ZHerm1_MR_Star );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A21, BHerm1_MR_Star, (T)1, ZHerm2_MC_Star );
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

        SlidePartitionDown
        ( ZHermT_MC_Star,  ZHerm0_MC_Star,
                           ZHerm1_MC_Star,
         /**************/ /**************/
          ZHermB_MC_Star,  ZHerm2_MC_Star );
        
        SlidePartitionDown
        ( ZHermT_MR_Star,  ZHerm0_MR_Star,
                           ZHerm1_MR_Star,
         /**************/ /**************/
          ZHermB_MR_Star,  ZHerm2_MR_Star );

    }

    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalHemmAccumulateRL
( T alpha,
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,Star,MC  >& B_Star_MC,
  const DistMatrix<T,Star,MR  >& B_Star_MR,
        DistMatrix<T,Star,MC  >& Z_Star_MC,
        DistMatrix<T,Star,MR  >& Z_Star_MR )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalHemmAccumulateRL");
    if( A.Grid() != B_Star_MC.Grid() ||
        B_Star_MC.Grid() != B_Star_MR.Grid() ||
        B_Star_MR.Grid() != Z_Star_MC.Grid() ||
        Z_Star_MC.Grid() != Z_Star_MR.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != A.Width() ||
        A.Height() != B_Star_MC.Width() ||
        A.Height() != B_Star_MR.Width() ||
        A.Height() != Z_Star_MC.Width() ||
        A.Height() != Z_Star_MR.Width() ||
        B_Star_MC.Height() != B_Star_MR.Height() ||
        B_Star_MR.Height() != Z_Star_MC.Height() ||
        Z_Star_MC.Height() != Z_Star_MR.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalHemmAccumulateRL: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B[* ,MC] ~ " << B_Star_MC.Height() << " x " 
                               << B_Star_MC.Width() << endl
            << "  B[* ,MR] ~ " << B_Star_MR.Height() << " x " 
                               << B_Star_MR.Width() << endl
            << "  Z[* ,MC] ~ " << Z_Star_MC.Height() << " x " 
                               << Z_Star_MC.Width() << endl
            << "  Z[* ,MR] ~ " << Z_Star_MR.Height() << " x " 
                               << Z_Star_MR.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( B_Star_MC.RowAlignment() != A.ColAlignment() ||
        B_Star_MR.RowAlignment() != A.RowAlignment() ||
        Z_Star_MC.RowAlignment() != A.ColAlignment() ||
        Z_Star_MR.RowAlignment() != A.RowAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,Star,MC>
        BL_Star_MC(g), BR_Star_MC(g),
        B0_Star_MC(g), B1_Star_MC(g), B2_Star_MC(g);

    DistMatrix<T,Star,MR>
        BL_Star_MR(g), BR_Star_MR(g),
        B0_Star_MR(g), B1_Star_MR(g), B2_Star_MR(g);

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
    LockedPartitionRight( B_Star_MR,  BL_Star_MR, BR_Star_MR, 0 );
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

        LockedRepartitionRight
        ( BL_Star_MR, /**/ BR_Star_MR,
          B0_Star_MR, /**/ B1_Star_MR, B2_Star_MR );

        RepartitionRight
        ( ZL_Star_MC, /**/ ZR_Star_MC,
          Z0_Star_MC, /**/ Z1_Star_MC, Z2_Star_MC );

        RepartitionRight
        ( ZL_Star_MR, /**/ ZR_Star_MR,
          Z0_Star_MR, /**/ Z1_Star_MR, Z2_Star_MR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Lower );
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, B1_Star_MC, D11, (T)1, Z1_Star_MR );
        D11.MakeTrapezoidal( Left, Lower, -1 );

        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose,
          alpha, B1_Star_MR, D11, (T)1, Z1_Star_MC );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, B2_Star_MC, A21, (T)1, Z1_Star_MR );

        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose,
          alpha, B1_Star_MR, A21, (T)1, Z2_Star_MC );
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

        SlideLockedPartitionRight
        ( BL_Star_MR,             /**/ BR_Star_MR,
          B0_Star_MR, B1_Star_MR, /**/ B2_Star_MR );

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


template void elemental::blas::internal::HemmRL
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::HemmRLA
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::LocalHemmAccumulateRL
( float alpha,
  const DistMatrix<float,MC,  MR  >& A,
  const DistMatrix<float,Star,MC  >& B_Star_MC,
  const DistMatrix<float,Star,MR  >& B_Star_MR,
        DistMatrix<float,Star,MC  >& Z_Star_MC,
        DistMatrix<float,Star,MR  >& Z_Star_MR );

template void elemental::blas::internal::HemmRL
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::HemmRLA
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::LocalHemmAccumulateRL
( double alpha,
  const DistMatrix<double,MC,  MR  >& A,
  const DistMatrix<double,Star,MC  >& B_Star_MC,
  const DistMatrix<double,Star,MR  >& B_Star_MR,
        DistMatrix<double,Star,MC  >& Z_Star_MC,
        DistMatrix<double,Star,MR  >& Z_Star_MR );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::HemmRL
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::HemmRLA
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::LocalHemmAccumulateRL
( scomplex alpha,
  const DistMatrix<scomplex,MC,  MR  >& A,
  const DistMatrix<scomplex,Star,MC  >& B_Star_MC,
  const DistMatrix<scomplex,Star,MR  >& B_Star_MR,
        DistMatrix<scomplex,Star,MC  >& Z_Star_MC,
        DistMatrix<scomplex,Star,MR  >& Z_Star_MR );

template void elemental::blas::internal::HemmRL
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::HemmRLA
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
    
template void elemental::blas::internal::LocalHemmAccumulateRL
( dcomplex alpha,
  const DistMatrix<dcomplex,MC,  MR  >& A,
  const DistMatrix<dcomplex,Star,MC  >& B_Star_MC,
  const DistMatrix<dcomplex,Star,MR  >& B_Star_MR,
        DistMatrix<dcomplex,Star,MC  >& Z_Star_MC,
        DistMatrix<dcomplex,Star,MR  >& Z_Star_MR );
#endif

