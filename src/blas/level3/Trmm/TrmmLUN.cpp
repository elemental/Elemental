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
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

#include "./TrmmUtil.hpp"

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

// Left Upper Normal (Non)Unit Trmm
//   X := triu(U)  X, or
//   X := triuu(U) X
template<typename T>
void
elemental::blas::internal::TrmmLUN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmLUN");
#endif
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Width() )
        blas::internal::TrmmLUNA( diagonal, alpha, U, X );
    else
        blas::internal::TrmmLUNC( diagonal, alpha, U, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::TrmmLUNA
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmULNA");
    if( U.Grid() != X.Grid() )
        throw logic_error( "U and X must be distributed over the same grid." );
    if( U.Height() != U.Width() || U.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLUNA: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = U.Grid();

    DistMatrix<T,MC,MR>
        XL(g), XR(g),
        X0(g), X1(g), X2(g);

    DistMatrix<T,VR,  Star> X1_VR_Star(g);
    DistMatrix<T,Star,MR  > X1Trans_Star_MR(g);
    DistMatrix<T,MC,  Star> Z1_MC_Star(g);

    PartitionRight
    ( X, XL, XR, 0 );
    while( XL.Width() < X.Width() )
    {
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        X1_VR_Star.AlignWith( U );
        X1Trans_Star_MR.AlignWith( U );
        Z1_MC_Star.AlignWith( U );
        Z1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VR_Star = X1;
        X1Trans_Star_MR.TransposeFrom( X1_VR_Star );
        Z1_MC_Star.SetToZero();
        blas::internal::LocalTrmmAccumulateLUN
        ( Transpose, diagonal, alpha, U, X1Trans_Star_MR, Z1_MC_Star );

        X1.SumScatterFrom( Z1_MC_Star );
        //--------------------------------------------------------------------//
        X1_VR_Star.FreeAlignments();
        X1Trans_Star_MR.FreeAlignments();
        Z1_MC_Star.FreeAlignments();

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::TrmmLUNC
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmLUNC");
    if( U.Grid() != X.Grid() )
        throw logic_error( "U and X must be distributed over the same grid." );
    if( U.Height() != U.Width() || U.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLUN: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T,MC,MR> XT(g),  X0(g),
                        XB(g),  X1(g),
                                X2(g);

    // Temporary distributions
    DistMatrix<T,Star,Star> U11_Star_Star(g);
    DistMatrix<T,Star,MC  > U12_Star_MC(g);
    DistMatrix<T,Star,VR  > X1_Star_VR(g);
    DistMatrix<T,Star,MR  > D1_Star_MR(g);

    // Start the algorithm
    blas::Scal( alpha, X );
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

        U12_Star_MC.AlignWith( X2 );
        D1_Star_MR.AlignWith( X1 );
        D1_Star_MR.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_Star_VR = X1;
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Left, Upper, Normal, diagonal, (T)1, U11_Star_Star, X1_Star_VR );
        X1 = X1_Star_VR;
 
        U12_Star_MC = U12;
        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, U12_Star_MC, X2, (T)0, D1_Star_MR );
        X1.SumScatterUpdate( (T)1, D1_Star_MR );
       //--------------------------------------------------------------------//
        U12_Star_MC.FreeAlignments();
        D1_Star_MR.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTrmmAccumulateLUN
( Orientation orientation, Diagonal diagonal, T alpha,
  const DistMatrix<T,MC,  MR  >& U,
  const DistMatrix<T,Star,MR  >& XHermOrTrans_Star_MR,
        DistMatrix<T,MC,  Star>& Z_MC_Star )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTrmmAccumulateLUN");
    if( U.Grid() != XHermOrTrans_Star_MR.Grid() ||
        XHermOrTrans_Star_MR.Grid() != Z_MC_Star.Grid() )
        throw logic_error( "{U,X,Z} must be distributed over the same grid." );
    if( U.Height() != U.Width() ||
        U.Height() != XHermOrTrans_Star_MR.Width() ||
        U.Height() != Z_MC_Star.Height() ||
        XHermOrTrans_Star_MR.Height() != Z_MC_Star.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateLUN: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X^H/T[* ,MR] ~ " << XHermOrTrans_Star_MR.Height() << " x "
                                   << XHermOrTrans_Star_MR.Width() << "\n"
            << "  Z[MC,* ] ~ " << Z_MC_Star.Height() << " x "
                               << Z_MC_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( XHermOrTrans_Star_MR.RowAlignment() != A.RowAlignment() ||
        Z_MC_Star.ColAlignment() != A.ColAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,Star,MR>
        XLHermOrTrans_Star_MR(g), XRHermOrTrans_Star_MR(g),
        X0HermOrTrans_Star_MR(g), X1HermOrTrans_Star_MR(g), 
        X2HermOrTrans_Star_MR(g);

    DistMatrix<T,MC,Star>
        ZT_MC_Star(g),  Z0_MC_Star(g),
        ZB_MC_Star(g),  Z1_MC_Star(g),
                        Z2_MC_Star(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    LockedPartitionRight
    ( XHermOrTrans_Star_MR, XLHermOrTrans_Star_MR, XRHermOrTrans_Star_MR, 0 );
    PartitionDown
    ( Z_MC_Star, ZT_MC_Star,
                 ZB_MC_Star, 0 );
    while( UTL.Height() < U.Height() )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        LockedRepartitionRight
        ( XLHermOrTrans_Star_MR, /**/ XRHermOrTrans_Star_MR,
          X0HermOrTrans_Star_MR, /**/ X1HermOrTrans_Star_MR, 
                                      X2HermOrTrans_Star_MR );

        RepartitionDown
        ( ZT_MC_Star,  Z0_MC_Star,
         /**********/ /**********/
                       Z1_MC_Star,
          ZB_MC_Star,  Z2_MC_Star );

        D11.AlignWith( U11 );
        //--------------------------------------------------------------------//
        D11 = U11;
        D11.MakeTrapezoidal( Left, Upper );
        if( diagonal == Unit )
            SetDiagonalToOne( D11 );
        blas::internal::LocalGemm
        ( Normal, orientation, alpha, D11, X1HermOrTrans_Star_MR,
          (T)1, Z1_MC_Star );

        blas::internal::LocalGemm
        ( Normal, orientation, alpha, U01, X1HermOrTrans_Star_MR,
          (T)1, Z0_MC_Star );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlideLockedPartitionRight
        ( XLHermOrTrans_Star_MR,                        /**/ 
          XRHermOrTrans_Star_MR,
          X0HermOrTrans_Star_MR, X1HermOrTrans_Star_MR, /**/ 
          X2HermOrTrans_Star_MR );

        SlidePartitionDown
        ( ZT_MC_Star,  Z0_MC_Star,
                       Z1_MC_Star,
         /**********/ /**********/
          ZB_MC_Star,  Z2_MC_Star );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrmmLUN
( Diagonal diagonal,
  float alpha, const DistMatrix<float,MC,MR>& U,
                     DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmLUN
( Diagonal diagonal,
  double alpha, const DistMatrix<double,MC,MR>& U,
                      DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmLUN
( Diagonal diagonal,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& U,
                        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmLUN
( Diagonal diagonal,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& U,
                        DistMatrix<dcomplex,MC,MR>& X );
#endif

