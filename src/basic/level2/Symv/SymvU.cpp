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
elemental::basic::internal::LocalSymvColAccumulateU
( T alpha, 
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,Star>& x_MC_Star,
  const DistMatrix<T,MR,Star>& x_MR_Star,
        DistMatrix<T,MC,Star>& z_MC_Star,
        DistMatrix<T,MR,Star>& z_MR_Star )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalSymvColAccumulateU");
    if( A.Grid() != x_MC_Star.Grid() ||
        x_MC_Star.Grid() != x_MR_Star.Grid() ||
        x_MR_Star.Grid() != z_MC_Star.Grid() ||
        z_MC_Star.Grid() != z_MR_Star.Grid() )
        throw logic_error( "{A,x,z} must be distributed over the same grid." );
    if( x_MC_Star.Width() != 1 || x_MR_Star.Width() != 1 ||
        z_MC_Star.Width() != 1 || z_MR_Star.Width() != 1 )
        throw logic_error( "Expected x and z to be column vectors." );
    if( A.Height() != A.Width() || 
        A.Height() != x_MC_Star.Height() ||
        A.Height() != x_MR_Star.Height() ||
        A.Height() != z_MC_Star.Height() ||
        A.Height() != z_MR_Star.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalSymvColAccumulateU: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x[MC,* ] ~ " << x_MC_Star.Height() << " x " 
                               << x_MC_Star.Width() << endl
            << "  x[MR,* ] ~ " << x_MR_Star.Height() << " x " 
                               << x_MR_Star.Width() << endl
            << "  z[MC,* ] ~ " << z_MC_Star.Height() << " x " 
                               << z_MC_Star.Width() << endl
            << "  z[MR,* ] ~ " << z_MR_Star.Height() << " x " 
                               << z_MR_Star.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( x_MC_Star.ColAlignment() != A.ColAlignment() ||
        x_MR_Star.ColAlignment() != A.RowAlignment() ||
        z_MC_Star.ColAlignment() != A.ColAlignment() ||
        z_MR_Star.ColAlignment() != A.RowAlignment() )
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
        xT_MC_Star(g),  x0_MC_Star(g),
        xB_MC_Star(g),  x1_MC_Star(g),
                        x2_MC_Star(g);

    DistMatrix<T,MR,Star> 
        xT_MR_Star(g),  x0_MR_Star(g),
        xB_MR_Star(g),  x1_MR_Star(g),
                        x2_MR_Star(g);

    DistMatrix<T,MC,Star> 
        zT_MC_Star(g),  z0_MC_Star(g),
        zB_MC_Star(g),  z1_MC_Star(g),
                        z2_MC_Star(g);

    DistMatrix<T,MR,Star> 
        zT_MR_Star(g),  z0_MR_Star(g),
        zB_MR_Star(g),  z1_MR_Star(g),
                        z2_MR_Star(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );
                 
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDown
    ( x_MC_Star, xT_MC_Star,
                 xB_MC_Star, 0 );
    LockedPartitionDown
    ( x_MR_Star, xT_MR_Star,
                 xB_MR_Star, 0 );
    PartitionDown
    ( z_MC_Star, zT_MC_Star,
                 zB_MC_Star, 0 );
    PartitionDown
    ( z_MR_Star, zT_MR_Star,
                 zB_MR_Star, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown
        ( xT_MC_Star,  x0_MC_Star,
         /**********/ /**********/
                       x1_MC_Star,
          xB_MC_Star,  x2_MC_Star );

        LockedRepartitionDown
        ( xT_MR_Star,  x0_MR_Star,
         /**********/ /**********/
                       x1_MR_Star,
          xB_MR_Star,  x2_MR_Star );

        RepartitionDown
        ( zT_MC_Star,  z0_MC_Star,
         /**********/ /**********/
                       z1_MC_Star,
          zB_MC_Star,  z2_MC_Star );

        RepartitionDown
        ( zT_MR_Star,  z0_MR_Star,
         /**********/ /**********/
                       z1_MR_Star,
          zB_MR_Star,  z2_MR_Star );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Upper );
        basic::Gemv
        ( Normal, 
          alpha, D11.LockedLocalMatrix(), 
                 x1_MR_Star.LockedLocalMatrix(),
          (T)1,  z1_MC_Star.LocalMatrix() );
        D11.MakeTrapezoidal( Left, Upper, 1 );

        basic::Gemv
        ( Transpose,
          alpha, D11.LockedLocalMatrix(),
                 x1_MC_Star.LockedLocalMatrix(),
          (T)1,  z1_MR_Star.LocalMatrix() );

        basic::Gemv
        ( Normal,
          alpha, A12.LockedLocalMatrix(),
                 x2_MR_Star.LockedLocalMatrix(),
          (T)1,  z1_MC_Star.LocalMatrix() );

        basic::Gemv
        ( Transpose,
          alpha, A12.LockedLocalMatrix(),
                 x1_MC_Star.LockedLocalMatrix(),
          (T)1,  z2_MR_Star.LocalMatrix() );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown
        ( xT_MC_Star,  x0_MC_Star,
                       x1_MC_Star,
         /**********/ /**********/
          xB_MC_Star,  x2_MC_Star );

        SlideLockedPartitionDown
        ( xT_MR_Star,  x0_MR_Star,
                       x1_MR_Star,
         /**********/ /**********/
          xB_MR_Star,  x2_MR_Star );
        
        SlidePartitionDown
        ( zT_MC_Star,  z0_MC_Star,
                       z1_MC_Star,
         /**********/ /**********/
          zB_MC_Star,  z2_MC_Star );

        SlidePartitionDown
        ( zT_MR_Star,  z0_MR_Star,
                       z1_MR_Star,
         /**********/ /**********/
          zB_MR_Star,  z2_MR_Star );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::basic::internal::LocalSymvRowAccumulateU
( T alpha, 
  const DistMatrix<T,MC,  MR>& A,
  const DistMatrix<T,Star,MC>& x_Star_MC,
  const DistMatrix<T,Star,MR>& x_Star_MR,
        DistMatrix<T,Star,MC>& z_Star_MC,
        DistMatrix<T,Star,MR>& z_Star_MR )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalSymvRowAccumulateU");
    if( A.Grid() != x_Star_MC.Grid() ||
        x_Star_MC.Grid() != x_Star_MR.Grid() ||
        x_Star_MR.Grid() != z_Star_MC.Grid() ||
        z_Star_MC.Grid() != z_Star_MR.Grid() )
        throw logic_error( "{A,x,z} must be distributed over the same grid." );
    if( x_Star_MC.Height() != 1 || x_Star_MR.Height() != 1 ||
        z_Star_MC.Height() != 1 || z_Star_MR.Height() != 1 )
        throw logic_error( "Expected x and z to be row vectors." );
    if( A.Height() != A.Width() || 
        A.Height() != x_Star_MC.Width() ||
        A.Height() != x_Star_MR.Width() ||
        A.Height() != z_Star_MC.Width() ||
        A.Height() != z_Star_MR.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal LocalSymvRowAccumulateU: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x[* ,MC] ~ " << x_Star_MC.Height() << " x " 
                               << x_Star_MC.Width() << endl
            << "  x[* ,MR] ~ " << x_Star_MR.Height() << " x " 
                               << x_Star_MR.Width() << endl
            << "  z[* ,MC] ~ " << z_Star_MC.Height() << " x " 
                               << z_Star_MC.Width() << endl
            << "  z[* ,MR] ~ " << z_Star_MR.Height() << " x " 
                               << z_Star_MR.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( x_Star_MC.RowAlignment() != A.ColAlignment() ||
        x_Star_MR.RowAlignment() != A.RowAlignment() ||
        z_Star_MC.RowAlignment() != A.ColAlignment() ||
        z_Star_MR.RowAlignment() != A.RowAlignment() )
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
        xL_Star_MC(g), xR_Star_MC(g),
        x0_Star_MC(g), x1_Star_MC(g), x2_Star_MC(g);

    DistMatrix<T,Star,MR> 
        xL_Star_MR(g), xR_Star_MR(g),
        x0_Star_MR(g), x1_Star_MR(g), x2_Star_MR(g);

    DistMatrix<T,Star,MC> 
        zL_Star_MC(g), zR_Star_MC(g),
        z0_Star_MC(g), z1_Star_MC(g), z2_Star_MC(g);

    DistMatrix<T,Star,MR> 
        zL_Star_MR(g), zR_Star_MR(g),
        z0_Star_MR(g), z1_Star_MR(g), z2_Star_MR(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );
                 
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( x_Star_MC,  xL_Star_MC, xR_Star_MC, 0 );
    LockedPartitionRight( x_Star_MR,  xL_Star_MR, xR_Star_MR, 0 );
    PartitionRight( z_Star_MC,  zL_Star_MC, zR_Star_MC, 0 );
    PartitionRight( z_Star_MR,  zL_Star_MR, zR_Star_MR, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight
        ( xL_Star_MC, /**/ xR_Star_MC, 
          x0_Star_MC, /**/ x1_Star_MC, x2_Star_MC );

        LockedRepartitionRight
        ( xL_Star_MR, /**/ xR_Star_MR, 
          x0_Star_MR, /**/ x1_Star_MR, x2_Star_MR );

        RepartitionRight
        ( zL_Star_MC, /**/ zR_Star_MC,
          z0_Star_MC, /**/ z1_Star_MC, z2_Star_MC );

        RepartitionRight
        ( zL_Star_MR, /**/ zR_Star_MR,
          z0_Star_MR, /**/ z1_Star_MR, z2_Star_MR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( Left, Upper );
        basic::Gemv
        ( Normal, 
          alpha, D11.LockedLocalMatrix(), 
                 x1_Star_MR.LockedLocalMatrix(),
          (T)1,  z1_Star_MC.LocalMatrix() );
        D11.MakeTrapezoidal( Left, Upper, 1 );

        basic::Gemv
        ( Transpose,
          alpha, D11.LockedLocalMatrix(),
                 x1_Star_MC.LockedLocalMatrix(),
          (T)1,  z1_Star_MR.LocalMatrix() );

        basic::Gemv
        ( Normal,
          alpha, A12.LockedLocalMatrix(),
                 x2_Star_MR.LockedLocalMatrix(),
          (T)1,  z1_Star_MC.LocalMatrix() );
        basic::Gemv
        ( Transpose,
          alpha, A12.LockedLocalMatrix(),
                 x1_Star_MC.LockedLocalMatrix(),
          (T)1,  z2_Star_MR.LocalMatrix() );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight
        ( xL_Star_MC,             /**/ xR_Star_MC,
          x0_Star_MC, x1_Star_MC, /**/ x2_Star_MC );
        
        SlideLockedPartitionRight
        ( xL_Star_MR,             /**/ xR_Star_MR,
          x0_Star_MR, x1_Star_MR, /**/ x2_Star_MR );

        SlidePartitionRight
        ( zL_Star_MC,             /**/ zR_Star_MC,
          z0_Star_MC, z1_Star_MC, /**/ z2_Star_MC );

        SlidePartitionRight
        ( zL_Star_MR,             /**/ zR_Star_MR,
          z0_Star_MR, z1_Star_MR, /**/ z2_Star_MR );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::internal::LocalSymvColAccumulateU
( float alpha, 
  const DistMatrix<float,MC,MR  >& A,
  const DistMatrix<float,MC,Star>& x_MC_Star,    
  const DistMatrix<float,MR,Star>& x_MR_Star,
        DistMatrix<float,MC,Star>& z_MC_Star,
        DistMatrix<float,MR,Star>& z_MR_Star );

template void elemental::basic::internal::LocalSymvRowAccumulateU
( float alpha, 
  const DistMatrix<float,MC,  MR>& A,
  const DistMatrix<float,Star,MC>& x_Star_MC,
  const DistMatrix<float,Star,MR>& x_Star_MR,
        DistMatrix<float,Star,MC>& z_Star_MC,
        DistMatrix<float,Star,MR>& z_Star_MR );

template void elemental::basic::internal::LocalSymvColAccumulateU
( double alpha, 
  const DistMatrix<double,MC,MR  >& A,
  const DistMatrix<double,MC,Star>& x_MC_Star,
  const DistMatrix<double,MR,Star>& x_MR_Star,
        DistMatrix<double,MC,Star>& z_MC_Star,
        DistMatrix<double,MR,Star>& z_MR_Star );

template void elemental::basic::internal::LocalSymvRowAccumulateU
( double alpha, 
  const DistMatrix<double,MC,  MR>& A,
  const DistMatrix<double,Star,MC>& x_Star_MC,
  const DistMatrix<double,Star,MR>& x_Star_MR,
        DistMatrix<double,Star,MC>& z_Star_MC,
        DistMatrix<double,Star,MR>& z_Star_MR );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::LocalSymvColAccumulateU
( scomplex alpha,
  const DistMatrix<scomplex,MC,MR  >& A,
  const DistMatrix<scomplex,MC,Star>& x_MC_Star,
  const DistMatrix<scomplex,MR,Star>& x_MR_Star,
        DistMatrix<scomplex,MC,Star>& z_MC_Star,
        DistMatrix<scomplex,MR,Star>& z_MR_Star );

template void elemental::basic::internal::LocalSymvRowAccumulateU
( scomplex alpha,
  const DistMatrix<scomplex,MC,  MR>& A,
  const DistMatrix<scomplex,Star,MC>& x_Star_MC,
  const DistMatrix<scomplex,Star,MR>& x_Star_MR,
        DistMatrix<scomplex,Star,MC>& z_Star_MC,
        DistMatrix<scomplex,Star,MR>& z_Star_MR );

template void elemental::basic::internal::LocalSymvColAccumulateU
( dcomplex alpha,
  const DistMatrix<dcomplex,MC,MR  >& A,
  const DistMatrix<dcomplex,MC,Star>& x_MC_Star,
  const DistMatrix<dcomplex,MR,Star>& x_MR_Star,
        DistMatrix<dcomplex,MC,Star>& z_MC_Star,
        DistMatrix<dcomplex,MR,Star>& z_MR_Star );

template void elemental::basic::internal::LocalSymvRowAccumulateU
( dcomplex alpha,
  const DistMatrix<dcomplex,MC,  MR>& A,
  const DistMatrix<dcomplex,Star,MC>& x_Star_MC,
  const DistMatrix<dcomplex,Star,MR>& x_Star_MR,
        DistMatrix<dcomplex,Star,MC>& z_Star_MC,
        DistMatrix<dcomplex,Star,MR>& z_Star_MR );
#endif

