/*
   Copyright (c) 2009-2011, Jack Poulson
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
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

template<typename T>
void
elemental::basic::internal::HemmRU
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmRU");
#endif
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Height() )
        basic::internal::HemmRUA( alpha, A, B, beta, C );
    else
        basic::internal::HemmRUC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::HemmRUA
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmRUA");
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

    DistMatrix<T,MR,  Star> B1Herm_MR_Star(g);
    DistMatrix<T,VC,  Star> B1Herm_VC_Star(g);
    DistMatrix<T,Star,MC  > B1_Star_MC(g);
    DistMatrix<T,MC,  Star> Z1Herm_MC_Star(g);
    DistMatrix<T,MR,  Star> Z1Herm_MR_Star(g);
    DistMatrix<T,MC,  MR  > Z1Herm(g);
    DistMatrix<T,MR,  MC  > Z1Herm_MR_MC(g);

    Matrix<T> Z1Local;

    basic::Scal( beta, C );
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
        Z1Herm_MC_Star.AlignWith( A );
        Z1Herm_MR_Star.AlignWith( A );
        Z1Herm_MR_MC.AlignWith( C1 );
        Z1Herm_MC_Star.ResizeTo( C1.Width(), C1.Height() );
        Z1Herm_MR_Star.ResizeTo( C1.Width(), C1.Height() );
        //--------------------------------------------------------------------//
        B1Herm_MR_Star.ConjugateTransposeFrom( B1 );
        B1Herm_VC_Star = B1Herm_MR_Star;
        B1_Star_MC.ConjugateTransposeFrom( B1Herm_VC_Star );
        Z1Herm_MC_Star.SetToZero();
        Z1Herm_MR_Star.SetToZero();
        basic::internal::LocalSymmetricAccumulateRU
        ( ConjugateTranspose, alpha, A, B1_Star_MC, B1Herm_MR_Star, 
          Z1Herm_MC_Star, Z1Herm_MR_Star );

        Z1Herm.SumScatterFrom( Z1Herm_MC_Star );
        Z1Herm_MR_MC = Z1Herm;
        Z1Herm_MR_MC.SumScatterUpdate( (T)1, Z1Herm_MR_Star );
        basic::ConjTrans( Z1Herm_MR_MC.LockedLocalMatrix(), Z1Local );
        basic::Axpy( (T)1, Z1Local, C1.LocalMatrix() );
        //--------------------------------------------------------------------//
        B1Herm_MR_Star.FreeAlignments();
        B1Herm_VC_Star.FreeAlignments();
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
elemental::basic::internal::HemmRUC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmRUC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw logic_error( "{A,B,C} must be distributed on the same grid." );
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
    DistMatrix<T,MC,Star> B1_MC_Star(g);
    DistMatrix<T,VR,  Star> AColPan_VR_Star(g);
    DistMatrix<T,Star,MR  > AColPanHerm_Star_MR(g);
    DistMatrix<T,MR,  Star> ARowPanHerm_MR_Star(g);

    // Start the algorithm
    basic::Scal( beta, C );
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

        basic::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          alpha, B1_MC_Star, ARowPanHerm_MR_Star, (T)1, CRight );

        basic::internal::LocalGemm
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
elemental::basic::internal::LocalSymmetricAccumulateRU
( Orientation orientation, T alpha,
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,Star,MC  >& B_Star_MC,
  const DistMatrix<T,MR,  Star>& BHermOrTrans_MR_Star,
        DistMatrix<T,MC,  Star>& ZHermOrTrans_MC_Star,
        DistMatrix<T,MR,  Star>& ZHermOrTrans_MR_Star )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalSymmetricAccumulateRU");
    if( A.Grid() != B_Star_MC.Grid() ||
        B_Star_MC.Grid() != BHermOrTrans_MR_Star.Grid() ||
        BHermOrTrans_MR_Star.Grid() != ZHermOrTrans_MC_Star.Grid() ||
        ZHermOrTrans_MC_Star.Grid() != ZHermOrTrans_MR_Star.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != A.Width() ||
        A.Height() != B_Star_MC.Width() ||
        A.Height() != BHermOrTrans_MR_Star.Height() ||
        A.Height() != ZHermOrTrans_MC_Star.Height() ||
        A.Height() != ZHermOrTrans_MR_Star.Height() ||
        B_Star_MC.Height() != BHermOrTrans_MR_Star.Width() ||
        BHermOrTrans_MR_Star.Width() != ZHermOrTrans_MC_Star.Width() ||
        ZHermOrTrans_MC_Star.Width() != ZHermOrTrans_MR_Star.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalSymmetricAccumulateRU: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B[* ,MC] ~ " << B_Star_MC.Height() << " x "
                               << B_Star_MC.Width() << endl
            << "  B^H/T[MR,* ] ~ " << BHermOrTrans_MR_Star.Height() << " x "
                                   << BHermOrTrans_MR_Star.Width() << endl
            << "  Z^H/T[MC,* ] ~ " << ZHermOrTrans_MC_Star.Height() << " x "
                                   << ZHermOrTrans_MC_Star.Width() << endl
            << "  Z^H/T[MR,* ] ~ " << ZHermOrTrans_MR_Star.Height() << " x "
                                   << ZHermOrTrans_MR_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( B_Star_MC.RowAlignment() != A.ColAlignment() ||
        BHermOrTrans_MR_Star.ColAlignment() != A.RowAlignment() ||
        ZHermOrTrans_MC_Star.ColAlignment() != A.ColAlignment() ||
        ZHermOrTrans_MR_Star.ColAlignment() != A.RowAlignment() )
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
        BTHermOrTrans_MR_Star(g),  B0HermOrTrans_MR_Star(g),
        BBHermOrTrans_MR_Star(g),  B1HermOrTrans_MR_Star(g),
                                   B2HermOrTrans_MR_Star(g);

    DistMatrix<T,MC,Star>
        ZTHermOrTrans_MC_Star(g),  Z0HermOrTrans_MC_Star(g),
        ZBHermOrTrans_MC_Star(g),  Z1HermOrTrans_MC_Star(g),
                                   Z2HermOrTrans_MC_Star(g);

    DistMatrix<T,MR,Star>
        ZBHermOrTrans_MR_Star(g),  Z0HermOrTrans_MR_Star(g),
        ZTHermOrTrans_MR_Star(g),  Z1HermOrTrans_MR_Star(g),
                                   Z2HermOrTrans_MR_Star(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( B_Star_MC,  BL_Star_MC, BR_Star_MC, 0 );
    LockedPartitionDown
    ( BHermOrTrans_MR_Star, BTHermOrTrans_MR_Star,
                            BBHermOrTrans_MR_Star, 0 );
    PartitionDown
    ( ZHermOrTrans_MC_Star, ZTHermOrTrans_MC_Star,
                            ZBHermOrTrans_MC_Star, 0 );
    PartitionDown
    ( ZHermOrTrans_MR_Star, ZTHermOrTrans_MR_Star,
                            ZBHermOrTrans_MR_Star, 0 );
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
        ( BTHermOrTrans_MR_Star,  B0HermOrTrans_MR_Star,
         /*********************/ /*********************/
                                  B1HermOrTrans_MR_Star,
          BBHermOrTrans_MR_Star,  B2HermOrTrans_MR_Star );

        RepartitionDown
        ( ZTHermOrTrans_MC_Star,  Z0HermOrTrans_MC_Star,
         /*********************/ /*********************/
                                  Z1HermOrTrans_MC_Star,
          ZBHermOrTrans_MC_Star,  Z2HermOrTrans_MC_Star );

        RepartitionDown
        ( ZTHermOrTrans_MR_Star,  Z0HermOrTrans_MR_Star,
         /*********************/ /*********************/
                                  Z1HermOrTrans_MR_Star,
          ZBHermOrTrans_MR_Star,  Z2HermOrTrans_MR_Star );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Upper );
        basic::internal::LocalGemm
        ( orientation, orientation,
          alpha, D11, B1_Star_MC, (T)1, Z1HermOrTrans_MR_Star );
        D11.MakeTrapezoidal( Left, Upper, 1 );

        basic::internal::LocalGemm
        ( Normal, Normal, alpha, D11, B1HermOrTrans_MR_Star, 
          (T)1, Z1HermOrTrans_MC_Star );

        basic::internal::LocalGemm
        ( orientation, orientation,
          alpha, A12, B1_Star_MC, (T)1, Z2HermOrTrans_MR_Star );

        basic::internal::LocalGemm
        ( Normal, Normal, alpha, A12, B2HermOrTrans_MR_Star, 
          (T)1, Z1HermOrTrans_MC_Star );
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
        ( BTHermOrTrans_MR_Star,  B0HermOrTrans_MR_Star,
                                  B1HermOrTrans_MR_Star,
         /*********************/ /*********************/
          BBHermOrTrans_MR_Star,  B2HermOrTrans_MR_Star );

        SlidePartitionDown
        ( ZTHermOrTrans_MC_Star,  Z0HermOrTrans_MC_Star,
                                  Z1HermOrTrans_MC_Star,
         /*********************/ /*********************/
          ZBHermOrTrans_MC_Star,  Z2HermOrTrans_MC_Star );       
        
        SlidePartitionDown
        ( ZTHermOrTrans_MR_Star,  Z0HermOrTrans_MR_Star,
                                  Z1HermOrTrans_MR_Star,
         /*********************/ /*********************/
          ZBHermOrTrans_MR_Star,  Z2HermOrTrans_MR_Star ); 
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::internal::HemmRU
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateRU
( Orientation orientation, float alpha,
  const DistMatrix<float,MC,  MR  >& A,
  const DistMatrix<float,Star,MC  >& B_Star_MC,
  const DistMatrix<float,MR,  Star>& BHermOrTrans_MR_Star,
        DistMatrix<float,MC,  Star>& ZHermOrTrans_MC_Star,
        DistMatrix<float,MR,  Star>& ZHermOrTrans_MR_Star );

template void elemental::basic::internal::HemmRU
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateRU
( Orientation orientation, double alpha,
  const DistMatrix<double,MC,  MR  >& A,
  const DistMatrix<double,Star,MC  >& B_Star_MC,
  const DistMatrix<double,MR,  Star>& BHermOrTrans_MR_Star,
        DistMatrix<double,MC,  Star>& ZHermOrTrans_MC_Star,
        DistMatrix<double,MR,  Star>& ZHermOrTrans_MR_Star );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::HemmRU
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateRU
( Orientation orientation, scomplex alpha,
  const DistMatrix<scomplex,MC,  MR  >& A,
  const DistMatrix<scomplex,Star,MC  >& B_Star_MC,
  const DistMatrix<scomplex,MR,  Star>& BHermOrTrans_MR_Star,
        DistMatrix<scomplex,MC,  Star>& ZHermOrTrans_MC_Star,
        DistMatrix<scomplex,MR,  Star>& ZHermOrTrans_MR_Star );

template void elemental::basic::internal::HemmRU
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateRU
( Orientation orientation, dcomplex alpha,
  const DistMatrix<dcomplex,MC,  MR  >& A,
  const DistMatrix<dcomplex,Star,MC  >& B_Star_MC,
  const DistMatrix<dcomplex,MR,  Star>& BHermOrTrans_MR_Star,
        DistMatrix<dcomplex,MC,  Star>& ZHermOrTrans_MC_Star,
        DistMatrix<dcomplex,MR,  Star>& ZHermOrTrans_MR_Star );
#endif

