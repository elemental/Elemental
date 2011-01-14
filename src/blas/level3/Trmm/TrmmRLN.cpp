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

// Right Lower Normal (Non)Unit Trmm
//   X := X tril(L), and
//   X := X trilu(L)
template<typename T>
void
elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRLN");
#endif
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        blas::internal::TrmmRLNA( diagonal, alpha, L, X );
    else
        blas::internal::TrmmRLNC( diagonal, alpha, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::TrmmRLNA
( Diagonal diagonal, 
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRLNA");
    if( L.Grid() != X.Grid() )
        throw logic_error( "{L,X} must be distributed over the same grid." );
#endif
    const Grid& g = L.Grid();

    DistMatrix<T,MC,MR>
        XT(g),  X0(g),
        XB(g),  X1(g),
                X2(g);

    DistMatrix<T,Star,VC  > X1_Star_VC(g);
    DistMatrix<T,Star,MC  > X1_Star_MC(g);
    DistMatrix<T,MR,  Star> Z1Trans_MR_Star(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

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

        X1_Star_VC.AlignWith( L );
        X1_Star_MC.AlignWith( L );
        Z1Trans_MR_Star.AlignWith( L );
        Z1Trans_MR_MC.AlignWith( X1 );
        Z1Trans_MR_Star.ResizeTo( X1.Width(), X1.Height() );
        //--------------------------------------------------------------------//
        X1_Star_VC = X1;
        X1_Star_MC = X1_Star_VC;
        Z1Trans_MR_Star.SetToZero();
        blas::internal::LocalTrmmAccumulateRLN
        ( Transpose, diagonal, alpha, L, X1_Star_MC, Z1Trans_MR_Star );

        Z1Trans_MR_MC.SumScatterFrom( Z1Trans_MR_Star );
        blas::Trans( Z1Trans_MR_MC.LocalMatrix(), X1.LocalMatrix() );
        //--------------------------------------------------------------------//
        X1_Star_VC.FreeAlignments();
        X1_Star_MC.FreeAlignments();
        Z1Trans_MR_Star.FreeAlignments();
        Z1Trans_MR_MC.FreeAlignments();

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
elemental::blas::internal::TrmmRLNC
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE

    PushCallStack("blas::internal::TrmmRLNC");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmRLNC: \n"
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
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,MR,  Star> L21_MR_Star(g);
    DistMatrix<T,VC,  Star> X1_VC_Star(g);
    DistMatrix<T,MC,  Star> D1_MC_Star(g);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );
 
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        L21_MR_Star.AlignWith( X2 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        L11_Star_Star = L11;
        blas::internal::LocalTrmm
        ( Right, Lower, Normal, diagonal, (T)1, L11_Star_Star, X1_VC_Star );
        X1 = X1_VC_Star;
 
        L21_MR_Star = L21;
        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, X2, L21_MR_Star, (T)0, D1_MC_Star );
        X1.SumScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        L21_MR_Star.FreeAlignments();
        D1_MC_Star.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

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
elemental::blas::internal::LocalTrmmAccumulateRLN
( Orientation orientation, Diagonal diagonal, T alpha,
  const DistMatrix<T,MC,  MR  >& L,
  const DistMatrix<T,Star,MC  >& X_Star_MC,
        DistMatrix<T,MR,  Star>& ZHermOrTrans_MR_Star )
{
#ifndef RELEASE
    PushCallStack("blas::internal::LocalTrmmAccumulateRLN");
    if( L.Grid() != X_Star_MC.Grid() || 
        X_Star_MC.Grid() != ZHermOrTrans_MR_Star.Grid() )
        throw logic_error( "{L,X,Z} must be distributed over the same grid." );
    if( L.Height() != L.Width() ||
        L.Height() != X_Star_MC.Width() ||
        L.Height() != ZHermOrTrans_MR_Star.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateRLN: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X[* ,MC] ~ " << X_Star_MC.Height() << " x "
                               << X_Star_MC.Width() << "\n"
            << "  Z^H/T[MR,* ] ~ " << ZHermOrTrans_MR_Star.Height() << " x "
                                   << ZHermOrTrans_MR_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( X_Star_MC.RowAlignment() != L.ColAlignment() ||
        ZHermOrTrans_MR_Star.ColAlignment() != L.RowAlignment() )
        throw logic_error( "Partial matrix distributions are misaligned." );
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,Star,MC>
        XL_Star_MC(g), XR_Star_MC(g),
        X0_Star_MC(g), X1_Star_MC(g), X2_Star_MC(g);

    DistMatrix<T,MR,Star>
        ZTHermOrTrans_MR_Star(g),  Z0HermOrTrans_MR_Star(g),
        ZBHermOrTrans_MR_Star(g),  Z1HermOrTrans_MR_Star(g),
                                   Z2HermOrTrans_MR_Star(g);

    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionRight( X_Star_MC,  XL_Star_MC, XR_Star_MC, 0 );
    PartitionDown
    ( ZHermOrTrans_MR_Star, ZTHermOrTrans_MR_Star,
                            ZBHermOrTrans_MR_Star, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionRight
        ( XL_Star_MC, /**/ XR_Star_MC,
          X0_Star_MC, /**/ X1_Star_MC, X2_Star_MC );

        RepartitionDown
        ( ZTHermOrTrans_MR_Star,  Z0HermOrTrans_MR_Star,
         /*********************/ /*********************/
                                  Z1HermOrTrans_MR_Star,
          ZBHermOrTrans_MR_Star,  Z2HermOrTrans_MR_Star );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        D11.MakeTrapezoidal( Left, Lower );
        if( diagonal == Unit )
            SetDiagonalToOne( D11 );
        blas::internal::LocalGemm
        ( orientation, orientation,
          alpha, D11, X1_Star_MC, (T)1, Z1HermOrTrans_MR_Star );

        blas::internal::LocalGemm
        ( orientation, orientation,
          alpha, L21, X2_Star_MC, (T)1, Z1HermOrTrans_MR_Star );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlideLockedPartitionRight
        ( XL_Star_MC,             /**/ XR_Star_MC,
          X0_Star_MC, X1_Star_MC, /**/ X2_Star_MC );

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

template void elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  float alpha, const DistMatrix<float,MC,MR>& L,
                     DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  double alpha, const DistMatrix<double,MC,MR>& L,
                      DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& L,
                        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& L,
                        DistMatrix<dcomplex,MC,MR>& X );
#endif

