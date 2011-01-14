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

// Left Lower Normal (Non)Unit Trmm 
//   X := tril(L)  X, or
//   X := trilu(L) X
template<typename T>
void
elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmLLN");
#endif
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        blas::internal::TrmmLLNA( diagonal, alpha, L, X );
    else
        blas::internal::TrmmLLNC( diagonal, alpha, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::TrmmLLNA
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmLLNA");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLLNA: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = L.Grid();

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

        X1_VR_Star.AlignWith( L );
        X1Trans_Star_MR.AlignWith( L );
        Z1_MC_Star.AlignWith( L );
        Z1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VR_Star = X1;
        X1Trans_Star_MR.TransposeFrom( X1_VR_Star );
        Z1_MC_Star.SetToZero();
        blas::internal::LocalTrmmAccumulateLLN
        ( Transpose, diagonal, alpha, L, X1Trans_Star_MR, Z1_MC_Star );

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
elemental::blas::internal::TrmmLLNC
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmLLNC");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLLNC: \n"
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

    DistMatrix<T,MC,MR> XT(g),  X0(g),
                        XB(g),  X1(g),
                                X2(g);

    // Temporary distributions
    DistMatrix<T,Star,MC  > L10_Star_MC(g);
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,Star,VR  > X1_Star_VR(g);
    DistMatrix<T,Star,MR  > D1_Star_MR(g);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        L10_Star_MC.AlignWith( X0 );
        D1_Star_MR.AlignWith( X1 );
        D1_Star_MR.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11;
        X1_Star_VR = X1;
        blas::internal::LocalTrmm
        ( Left, Lower, Normal, diagonal, (T)1, L11_Star_Star, X1_Star_VR );
        X1 = X1_Star_VR;

        L10_Star_MC = L10;
        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, L10_Star_MC, X0, (T)0, D1_Star_MR );
        X1.SumScatterUpdate( (T)1, D1_Star_MR );
        //--------------------------------------------------------------------//
        L10_Star_MC.FreeAlignments();
        D1_Star_MR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12, 
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::LocalTrmmAccumulateLLN
( Orientation orientation, Diagonal diagonal, T alpha,
  const DistMatrix<T,MC,  MR  >& L,
  const DistMatrix<T,Star,MR  >& XHermOrTrans_Star_MR,
        DistMatrix<T,MC,  Star>& Z_MC_Star )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTrmmAccumulateLLN");
    if( L.Grid() != XHermOrTrans_Star_MR.Grid() ||
        XHermOrTrans_Star_MR.Grid() != Z_MC_Star.Grid() )
        throw logic_error( "{L,X,Z} must be distributed over the same grid." );
    if( L.Height() != L.Width() ||
        L.Height() != XHermOrTrans_Star_MR.Width() ||
        L.Height() != Z_MC_Star.Height() ||
        XHermOrTrans_Star_MR.Height() != Z_MC_Star.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateLLN: \n" 
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X^H/T[* ,MR] ~ " << XHermOrTrans_Star_MR.Height() << " x "
                                   << XHermOrTrans_Star_MR.Width() << "\n"
            << "  Z[MC,* ] ~ " << Z_MC_Star.Height() << " x "
                               << Z_MC_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( XHermOrTrans_Star_MR.RowAlignment() != L.RowAlignment() ||
        Z_MC_Star.ColAlignment() != L.ColAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

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
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionRight
    ( XHermOrTrans_Star_MR, XLHermOrTrans_Star_MR, XRHermOrTrans_Star_MR, 0 );
    PartitionDown
    ( Z_MC_Star, ZT_MC_Star,
                 ZB_MC_Star, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionRight
        ( XLHermOrTrans_Star_MR, /**/ XRHermOrTrans_Star_MR,
          X0HermOrTrans_Star_MR, /**/ X1HermOrTrans_Star_MR, 
                                      X2HermOrTrans_Star_MR );

        RepartitionDown
        ( ZT_MC_Star,  Z0_MC_Star,
         /**********/ /**********/
                       Z1_MC_Star,
          ZB_MC_Star,  Z2_MC_Star );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        D11.MakeTrapezoidal( Left, Lower );
        if( diagonal == Unit )
            SetDiagonalToOne( D11 );
        blas::internal::LocalGemm
        ( Normal, orientation, alpha, D11, X1HermOrTrans_Star_MR,
          (T)1, Z1_MC_Star );

        blas::internal::LocalGemm
        ( Normal, orientation, alpha, L21, X1HermOrTrans_Star_MR,
          (T)1, Z2_MC_Star );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

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

template void elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  float alpha, const DistMatrix<float,MC,MR>& L,
                     DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  double alpha, const DistMatrix<double,MC,MR>& L,
                      DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& L,
                        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& L,
                        DistMatrix<dcomplex,MC,MR>& X );
#endif

