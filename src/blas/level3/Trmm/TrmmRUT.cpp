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

// Right Upper (Conjugate)Transpose (Non)Unit Trmm
//   X := X triu(U)^T, 
//   X := X triu(U)^H,
//   X := X triuu(U)^T, or
//   X := X triuu(U)^H
template<typename T>
void
elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRUT");
#endif
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        blas::internal::TrmmRUTA( orientation, diagonal, alpha, U, X );
    else
        blas::internal::TrmmRUTC( orientation, diagonal, alpha, U, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::TrmmRUTA
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRUTA");
    if( U.Grid() != X.Grid() )
        throw logic_error( "{U,X} must be distributed over the same grid." );
#endif
    const Grid& g = U.Grid();

    DistMatrix<T,MC,MR>
        XT(g),  X0(g),
        XB(g),  X1(g),
                X2(g);

    DistMatrix<T,MR,  Star> X1HermOrTrans_MR_Star(g);
    DistMatrix<T,MC,  Star> Z1HermOrTrans_MC_Star(g);
    DistMatrix<T,MC,  MR  > Z1HermOrTrans(g);
    DistMatrix<T,MR,  MC  > Z1HermOrTrans_MR_MC(g);

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

        X1HermOrTrans_MR_Star.AlignWith( U );
        Z1HermOrTrans_MC_Star.AlignWith( U );
        Z1HermOrTrans_MR_MC.AlignWith( X1 );
        Z1HermOrTrans_MC_Star.ResizeTo( X1.Width(), X1.Height() );
        //--------------------------------------------------------------------//
        if( orientation == ConjugateTranspose )
            X1HermOrTrans_MR_Star.ConjugateTransposeFrom( X1 );
        else
            X1HermOrTrans_MR_Star.TransposeFrom( X1 );
        Z1HermOrTrans_MC_Star.SetToZero();
        blas::internal::LocalTrmmAccumulateRUT
        ( orientation, diagonal, alpha,
          U, X1HermOrTrans_MR_Star, Z1HermOrTrans_MC_Star );

        Z1HermOrTrans.SumScatterFrom( Z1HermOrTrans_MC_Star );
        Z1HermOrTrans_MR_MC = Z1HermOrTrans;
        if( orientation == ConjugateTranspose )
        {
            blas::ConjTrans
            ( Z1HermOrTrans_MR_MC.LocalMatrix(), X1.LocalMatrix() );
        }
        else
        {
            blas::Trans
            ( Z1HermOrTrans_MR_MC.LocalMatrix(), X1.LocalMatrix() );
        }
        //--------------------------------------------------------------------//
        X1HermOrTrans_MR_Star.FreeAlignments();
        Z1HermOrTrans_MC_Star.FreeAlignments();
        Z1HermOrTrans_MR_MC.FreeAlignments();

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
elemental::blas::internal::TrmmRUTC
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRUTC");
    if( U.Grid() != X.Grid() )
        throw logic_error( "U and X must be distributed over the same grid." );
    if( orientation == Normal )
        throw logic_error( "TrmmRUTC expects a (Conjugate)Transpose option." );
    if( U.Height() != U.Width() || X.Width() != U.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmRUTC: \n"
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

    DistMatrix<T,MC,MR> XL(g), XR(g),
                        X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,Star,Star> U11_Star_Star(g);
    DistMatrix<T,Star,MR  > U12_Star_MR(g);
    DistMatrix<T,VC,  Star> X1_VC_Star(g);
    DistMatrix<T,MC,  Star> D1_MC_Star(g);
    
    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        U12_Star_MR.AlignWith( X2 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Right, Upper, orientation, diagonal, 
          (T)1, U11_Star_Star, X1_VC_Star );
        X1 = X1_VC_Star;
 
        U12_Star_MR = U12;
        blas::internal::LocalGemm
        ( Normal, orientation, (T)1, X2, U12_Star_MR, (T)0, D1_MC_Star );
        X1.SumScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        U12_Star_MR.FreeAlignments();
        D1_MC_Star.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12, 
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

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
elemental::blas::internal::LocalTrmmAccumulateRUT
( Orientation orientation, Diagonal diagonal, T alpha,
  const DistMatrix<T,MC,MR  >& U,
  const DistMatrix<T,MR,Star>& XHermOrTrans_MR_Star,
        DistMatrix<T,MC,Star>& ZHermOrTrans_MC_Star )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTrmmAccumulateRUT");
    if( U.Grid() != XHermOrTrans_MR_Star.Grid() ||
        XHermOrTrans_MR_Star.Grid() != ZHermOrTrans_MC_Star.Grid() )
        throw logic_error( "{U,X,Z} must be distributed over the same grid." );
    if( U.Height() != U.Width() ||
        U.Height() != XHermOrTrans_MR_Star.Height() ||
        U.Height() != ZHermOrTrans_MC_Star.Height() ||
        XHermOrTrans_MR_Star.Width() != ZHermOrTrans_MC_Star.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateRUT: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X^H/T[MR,* ] ~ " << XHermOrTrans_MR_Star.Height() << " x "
                                   << XHermOrTrans_MR_Star.Width() << "\n"
            << "  Z^H/T[MC,* ] ~ " << ZHermOrTrans_MC_Star.Height() << " x "
                                   << ZHermOrTrans_MC_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( XHermOrTrans_MR_Star.ColAlignment() != U.RowAlignment() ||
        ZHermOrTrans_MC_Star.ColAlignment() != U.ColAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,MR,Star>
        XTHermOrTrans_MR_Star(g),  X0HermOrTrans_MR_Star(g),
        XBHermOrTrans_MR_Star(g),  X1HermOrTrans_MR_Star(g),
                                   X2HermOrTrans_MR_Star(g);

    DistMatrix<T,MC,Star>
        ZTHermOrTrans_MC_Star(g),  Z0HermOrTrans_MC_Star(g),
        ZBHermOrTrans_MC_Star(g),  Z1HermOrTrans_MC_Star(g),
                                   Z2HermOrTrans_MC_Star(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    LockedPartitionDown
    ( XHermOrTrans_MR_Star, XTHermOrTrans_MR_Star,
                            XBHermOrTrans_MR_Star, 0 );
    PartitionDown
    ( ZHermOrTrans_MC_Star, ZTHermOrTrans_MC_Star,
                            ZBHermOrTrans_MC_Star, 0 );
    while( UTL.Height() < U.Height() )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        LockedRepartitionDown
        ( XTHermOrTrans_MR_Star,  X0HermOrTrans_MR_Star,
         /*********************/ /*********************/
                                  X1HermOrTrans_MR_Star,
          XBHermOrTrans_MR_Star,  X2HermOrTrans_MR_Star );

        RepartitionDown
        ( ZTHermOrTrans_MC_Star,  Z0HermOrTrans_MC_Star,
         /*********************/ /*********************/
                                  Z1HermOrTrans_MC_Star,
          ZBHermOrTrans_MC_Star,  Z2HermOrTrans_MC_Star );

        D11.AlignWith( U11 );
        //--------------------------------------------------------------------//
        D11 = U11;
        D11.MakeTrapezoidal( Left, Upper );
        if( diagonal == Unit )
            SetDiagonalToOne( D11 );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, D11, X1HermOrTrans_MR_Star,
          (T)1, Z1HermOrTrans_MC_Star );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, U01, X1HermOrTrans_MR_Star,
          (T)1, Z0HermOrTrans_MC_Star );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlideLockedPartitionDown
        ( XTHermOrTrans_MR_Star,  X0HermOrTrans_MR_Star,
                                  X1HermOrTrans_MR_Star,
         /*********************/ /*********************/
          XBHermOrTrans_MR_Star,  X2HermOrTrans_MR_Star );

        SlidePartitionDown
        ( ZTHermOrTrans_MC_Star,  Z0HermOrTrans_MC_Star,
                                  Z1HermOrTrans_MC_Star,
         /*********************/ /*********************/
          ZBHermOrTrans_MC_Star,  Z2HermOrTrans_MC_Star );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

