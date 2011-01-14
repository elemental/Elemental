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

// Right Lower (Conjugate)Transpose (Non)Unit Trmm
//   X := X tril(L)^T,
//   X := X tril(L)^H,
//   X := X trilu(L)^T, or
//   X := X trilu(L)^H
template<typename T>
void
elemental::blas::internal::TrmmRLT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRLT");
#endif
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        blas::internal::TrmmRLTA( orientation, diagonal, alpha, L, X );
    else
        blas::internal::TrmmRLTC( orientation, diagonal, alpha, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::TrmmRLTA
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRLTA");
    if( L.Grid() != X.Grid() )
        throw logic_error( "{L,X} must be distributed over the same grid." );
#endif
    const Grid& g = L.Grid();

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

        X1HermOrTrans_MR_Star.AlignWith( L );
        Z1HermOrTrans_MC_Star.AlignWith( L );
        Z1HermOrTrans_MR_MC.AlignWith( X1 );
        Z1HermOrTrans_MC_Star.ResizeTo( X1.Width(), X1.Height() );
        //--------------------------------------------------------------------//
        if( orientation == ConjugateTranspose )
            X1HermOrTrans_MR_Star.ConjugateTransposeFrom( X1 );
        else
            X1HermOrTrans_MR_Star.TransposeFrom( X1 );
        Z1HermOrTrans_MC_Star.SetToZero();
        blas::internal::LocalTrmmAccumulateRLT
        ( orientation, diagonal, alpha, 
          L, X1HermOrTrans_MR_Star, Z1HermOrTrans_MC_Star );

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
elemental::blas::internal::TrmmRLTC
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRLTC");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( orientation == Normal )
        throw logic_error( "TrmmRLTC expects a (Conjugate)Transpose option." );
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmRLTC: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,MC,MR> XL(g), XR(g),
                        X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,Star,MR  > L10_Star_MR(g);
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,VC,  Star> X1_VC_Star(g);
    DistMatrix<T,MC,  Star> D1_MC_Star(g);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionLeft( X, XL, XR, 0 );
    while( XL.Width() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionLeft
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );

        L10_Star_MR.AlignWith( X0 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        L11_Star_Star = L11;
        blas::internal::LocalTrmm
        ( Right, Lower, orientation, diagonal, 
          (T)1, L11_Star_Star, X1_VC_Star );
        X1 = X1_VC_Star;
 
        L10_Star_MR = L10;
        blas::internal::LocalGemm
        ( Normal, orientation, (T)1, X0, L10_Star_MR, (T)0, D1_MC_Star );
        X1.SumScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        L10_Star_MR.FreeAlignments();
        D1_MC_Star.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionLeft
        ( XL, /**/     XR,
          X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTrmmAccumulateRLT
( Orientation orientation, Diagonal diagonal, T alpha,
  const DistMatrix<T,MC,MR  >& L,
  const DistMatrix<T,MR,Star>& XHermOrTrans_MR_Star,
        DistMatrix<T,MC,Star>& ZHermOrTrans_MC_Star )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTrmmAccumulateRLT");
    if( L.Grid() != XHermOrTrans_MR_Star.Grid() ||
        XHermOrTrans_MR_Star.Grid() != ZHermOrTrans_MC_Star.Grid() )
        throw logic_error( "{L,X,Z} must be distributed over the same grid." );
    if( L.Height() != L.Width() ||
        L.Height() != XHermOrTrans_MR_Star.Height() ||
        L.Height() != ZHermOrTrans_MC_Star.Height() ||
        BHermOrTrans_MR_Star.Width() != ZHermOrTrans_MC_Star.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateRLT: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X^H/T[MR,* ] ~ " << XHermOrTrans_MR_Star.Height() << " x "
                                   << XHermOrTrans_MR_Star.Width() << "\n"
            << "  Z^H/T[MC,* ] ~ " << ZHermOrTrans_MC_Star.Height() << " x "
                                   << ZHermOrTrans_MC_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( XHermOrTrans_MR_Star.ColAlignment() != L.RowAlignment() ||
        ZHermOrTrans_MC_Star.ColAlignment() != L.ColAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

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
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionDown
    ( XHermOrTrans_MR_Star, XTHermOrTrans_MR_Star,
                            XBHermOrTrans_MR_Star, 0 );
    PartitionDown
    ( ZHermOrTrans_MC_Star, ZTHermOrTrans_MC_Star,
                            ZBHermOrTrans_MC_Star, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

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

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        D11.MakeTrapezoidal( Left, Lower );
        if( diagonal == Unit )
            SetDiagonalToOne( D11 );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, D11, X1HermOrTrans_MR_Star, 
          (T)1, Z1HermOrTrans_MC_Star );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, L21, X1HermOrTrans_MR_Star, 
          (T)1, Z2HermOrTrans_MC_Star );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

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

template void elemental::blas::internal::TrmmRLT
( Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmRLT
( Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmRLT
( Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmRLT
( Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

