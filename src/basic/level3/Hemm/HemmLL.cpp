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
elemental::basic::internal::HemmLL
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmLL");
#endif
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Width() )
        basic::internal::HemmLLA( alpha, A, B, beta, C );
    else
        basic::internal::HemmLLC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::HemmLLA
( T alpha, const DistMatrix<T,MC,MR>& A, 
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmLLA");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.Grid();

    DistMatrix<T,MC,MR> 
        BL(g), BR(g),
        B0(g), B1(g), B2(g);

    DistMatrix<T,MC,MR>
        CL(g), CR(g),
        C0(g), C1(g), C2(g);

    DistMatrix<T,MC,Star> B1_MC_Star(g);
    DistMatrix<T,VR,Star> B1_VR_Star(g);
    DistMatrix<T,Star,MR> B1Herm_Star_MR(g);
    DistMatrix<T,MC,MR  > Z1(g);
    DistMatrix<T,MC,Star> Z1_MC_Star(g);
    DistMatrix<T,MR,Star> Z1_MR_Star(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    basic::Scal( beta, C );
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
        B1_VR_Star.AlignWith( A );
        B1Herm_Star_MR.AlignWith( A );
        Z1_MC_Star.AlignWith( A );
        Z1_MR_Star.AlignWith( A );
        Z1.AlignWith( C1 );
        Z1_MC_Star.ResizeTo( C1.Height(), C1.Width() );
        Z1_MR_Star.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        B1_MC_Star = B1;
        B1_VR_Star = B1_MC_Star;
        B1Herm_Star_MR.ConjugateTransposeFrom( B1_VR_Star );
        Z1_MC_Star.SetToZero();
        Z1_MR_Star.SetToZero();
        basic::internal::LocalSymmetricAccumulateLL
        ( ConjugateTranspose, 
          alpha, A, B1_MC_Star, B1Herm_Star_MR, Z1_MC_Star, Z1_MR_Star );

        Z1_MR_MC.SumScatterFrom( Z1_MR_Star );
        Z1 = Z1_MR_MC;
        Z1.SumScatterUpdate( (T)1, Z1_MC_Star );
        basic::Axpy( (T)1, Z1, C1 );
        //--------------------------------------------------------------------//
        B1_MC_Star.FreeAlignments();
        B1_VR_Star.FreeAlignments();
        B1Herm_Star_MR.FreeAlignments();
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
elemental::basic::internal::HemmLLC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HemmLLC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.Grid();

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
    DistMatrix<T,MC,  Star> AColPan_MC_Star(g);
    DistMatrix<T,Star,MC  > ARowPan_Star_MC(g);
    DistMatrix<T,MR,  Star> B1Herm_MR_Star(g);

    // Start the algorithm
    basic::Scal( beta, C );
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
        B1Herm_MR_Star.AlignWith( C );
        //--------------------------------------------------------------------//
        AColPan_MC_Star = AColPan;
        ARowPan_Star_MC = ARowPan;
        AColPan_MC_Star.MakeTrapezoidal( Left, Lower );
        ARowPan_Star_MC.MakeTrapezoidal( Right, Lower, -1 );

        B1Herm_MR_Star.ConjugateTransposeFrom( B1 );

        basic::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          alpha, AColPan_MC_Star, B1Herm_MR_Star, (T)1, CBelow );

        basic::internal::LocalGemm
        ( ConjugateTranspose, ConjugateTranspose, 
          alpha, ARowPan_Star_MC, B1Herm_MR_Star, (T)1, CAbove );
        //--------------------------------------------------------------------//
        AColPan_MC_Star.FreeAlignments();
        ARowPan_Star_MC.FreeAlignments();
        B1Herm_MR_Star.FreeAlignments();

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
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, T alpha,  
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,MC,  Star>& B_MC_Star,
  const DistMatrix<T,Star,MR  >& BHermOrTrans_Star_MR,
        DistMatrix<T,MC,  Star>& Z_MC_Star,
        DistMatrix<T,MR,  Star>& Z_MR_Star )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalSymmetricAccumulateLL");
    if( A.Grid() != B_MC_Star.Grid() || 
        B_MC_Star.Grid() != BHermOrTrans_Star_MR.Grid() ||
        BHermOrTrans_Star_MR.Grid() != Z_MC_Star.Grid() ||
        Z_MC_Star.Grid() != Z_MR_Star.Grid() )
        throw logic_error( "{A,B,Z} must be distributed over the same grid." );
    if( A.Height() != A.Width() || 
        A.Height() != B_MC_Star.Height() ||
        A.Height() != BHermOrTrans_Star_MR.Width() ||
        A.Height() != Z_MC_Star.Height() ||
        A.Height() != Z_MR_Star.Height() ||
        B_MC_Star.Width() != BHermOrTrans_Star_MR.Height() ||
        BHermOrTrans_Star_MR.Height() != Z_MC_Star.Width() ||
        Z_MC_Star.Width() != Z_MR_Star.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalSymmetricAccumulateLL: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B[MC,* ] ~ " << B_MC_Star.Height() << " x " 
                               << B_MC_Star.Width() << "\n"
            << "  B^H[* ,MR] ~ " << BHermOrTrans_Star_MR.Height() << " x " 
                               << BHermOrTrans_Star_MR.Width() << "\n"
            << "  Z[MC,* ] ~ " << Z_MC_Star.Height() << " x " 
                               << Z_MC_Star.Width() << "\n"
            << "  Z[MR,* ] ` " << Z_MR_Star.Height() << " x " 
                               << Z_MR_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( B_MC_Star.ColAlignment() != A.ColAlignment() || 
        BHermOrTrans_Star_MR.RowAlignment() != A.RowAlignment() ||
        Z_MC_Star.ColAlignment() != A.ColAlignment() ||
        Z_MR_Star.ColAlignment() != A.RowAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = A.Grid();

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

    DistMatrix<T,Star,MR>
        BLHermOrTrans_Star_MR(g), BRHermOrTrans_Star_MR(g),
        B0HermOrTrans_Star_MR(g), B1HermOrTrans_Star_MR(g), 
        B2HermOrTrans_Star_MR(g);

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
    LockedPartitionRight
    ( BHermOrTrans_Star_MR, BLHermOrTrans_Star_MR, BRHermOrTrans_Star_MR, 0 );
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

        LockedRepartitionRight
        ( BLHermOrTrans_Star_MR, /**/ BRHermOrTrans_Star_MR,
          B0HermOrTrans_Star_MR, /**/ B1HermOrTrans_Star_MR, 
                                      B2HermOrTrans_Star_MR );

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
        basic::internal::LocalGemm
        ( Normal, orientation, alpha, D11, B1HermOrTrans_Star_MR, 
          (T)1, Z1_MC_Star );
        D11.MakeTrapezoidal( Left, Lower, -1 );

        basic::internal::LocalGemm
        ( orientation, Normal, 
          alpha, D11, B1_MC_Star, (T)1, Z1_MR_Star );

        basic::internal::LocalGemm
        ( Normal, orientation, alpha, A21, B1HermOrTrans_Star_MR, 
          (T)1, Z2_MC_Star );

        basic::internal::LocalGemm
        ( orientation, Normal, 
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

        SlideLockedPartitionRight
        ( BLHermOrTrans_Star_MR,                        /**/ 
          BRHermOrTrans_Star_MR,
          B0HermOrTrans_Star_MR, B1HermOrTrans_Star_MR, /**/ 
          B2HermOrTrans_Star_MR );

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

template void elemental::basic::internal::HemmLL
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, float alpha,  
  const DistMatrix<float,MC,  MR  >& A,
  const DistMatrix<float,MC,  Star>& B_MC_Star,
  const DistMatrix<float,Star,MR  >& BHermOrTrans_Star_MR,
        DistMatrix<float,MC,  Star>& Z_MC_Star,
        DistMatrix<float,MR,  Star>& Z_MR_Star );

template void elemental::basic::internal::HemmLL
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, double alpha,  
  const DistMatrix<double,MC,  MR  >& A,
  const DistMatrix<double,MC,  Star>& B_MC_Star,
  const DistMatrix<double,Star,MR  >& BHermOrTrans_Star_MR,
        DistMatrix<double,MC,  Star>& Z_MC_Star,
        DistMatrix<double,MR,  Star>& Z_MR_Star );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::HemmLL
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, scomplex alpha,  
  const DistMatrix<scomplex,MC,  MR  >& A,
  const DistMatrix<scomplex,MC,  Star>& B_MC_Star,
  const DistMatrix<scomplex,Star,MR  >& BHermOrTrans_Star_MR,
        DistMatrix<scomplex,MC,  Star>& Z_MC_Star,
        DistMatrix<scomplex,MR,  Star>& Z_MR_Star );

template void elemental::basic::internal::HemmLL
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void
elemental::basic::internal::LocalSymmetricAccumulateLL
( Orientation orientation, dcomplex alpha,  
  const DistMatrix<dcomplex,MC,  MR  >& A,
  const DistMatrix<dcomplex,MC,  Star>& B_MC_Star,
  const DistMatrix<dcomplex,Star,MR  >& BHermOrTrans_Star_MR,
        DistMatrix<dcomplex,MC,  Star>& Z_MC_Star,
        DistMatrix<dcomplex,MR,  Star>& Z_MR_Star );
#endif

