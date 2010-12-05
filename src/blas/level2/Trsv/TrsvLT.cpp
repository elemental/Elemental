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
elemental::blas::internal::TrsvLT
( Orientation orientation,
  Diagonal diagonal, 
  const DistMatrix<T,MC,MR>& L, 
        DistMatrix<T,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrsvLT");
    if( L.Grid() != x.Grid() )
        throw logic_error( "L and x must be distributed over the same grid." );
    if( orientation == Normal )
        throw logic_error( "TrsvLT expects a (conjugate-)transpose option." );
    if( L.Height() != L.Width() )
        throw logic_error( "L must be square." );
    if( x.Width() != 1 && x.Height() != 1 )
        throw logic_error( "x must be a vector." );
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( L.Width() != xLength )
        throw logic_error( "Nonconformal TrsvLT." );
#endif
    const Grid& g = L.Grid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<T,MC,MR> 
            LTL(g), LTR(g),  L00(g), L01(g), L02(g),
            LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                             L20(g), L21(g), L22(g);

        DistMatrix<T,MC,MR> 
            xT(g),  x0(g),
            xB(g),  x1(g),
                    x2(g);

        // Temporary distributions
        DistMatrix<T,Star,Star> L11_Star_Star(g);
        DistMatrix<T,Star,Star> x1_Star_Star(g);
        DistMatrix<T,MC,  Star> x1_MC_Star(g);
        DistMatrix<T,MR,  Star> z0_MR_Star(g);
        DistMatrix<T,MR,  MC  > z0_MR_MC(g);
        DistMatrix<T,MC,  MR  > z0(g);

        // Start the algorithm
        LockedPartitionUpDiagonal
        ( L, LTL, LTR,
             LBL, LBR, 0 );
        PartitionUp
        ( x, xT,
             xB, 0 );
        while( xT.Height() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            RepartitionUp
            ( xT,  x0,
                   x1,
             /**/ /**/
              xB,  x2 );

            x1_MC_Star.AlignWith( L10 );
            z0_MR_Star.AlignWith( L10 );
            z0_MR_Star.ResizeTo( x0.Height(), 1 );
            z0.AlignWith( x0 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            blas::Trsv
            ( Lower, orientation, diagonal,
              L11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_MC_Star = x1_Star_Star;
            blas::Gemv
            ( orientation, (T)-1, 
              L10.LockedLocalMatrix(), 
              x1_MC_Star.LockedLocalMatrix(),
              (T)0, z0_MR_Star.LocalMatrix() );
            z0_MR_MC.SumScatterFrom( z0_MR_Star );
            z0 = z0_MR_MC;
            blas::Axpy( (T)1, z0, x0 );
            //----------------------------------------------------------------//
            x1_MC_Star.FreeAlignments();
            z0_MR_Star.FreeAlignments();
            z0.FreeAlignments();

            SlideLockedPartitionUpDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            SlidePartitionUp
            ( xT,  x0,
             /**/ /**/
                   x1,
              xB,  x2 );
        }
    }
    else
    {
        // Matrix views 
        DistMatrix<T,MC,MR> 
            LTL(g), LTR(g),  L00(g), L01(g), L02(g),
            LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                             L20(g), L21(g), L22(g);

        DistMatrix<T,MC,MR> 
            xL(g), xR(g),
            x0(g), x1(g), x2(g);

        // Temporary distributions
        DistMatrix<T,Star,Star> L11_Star_Star(g);
        DistMatrix<T,Star,Star> x1_Star_Star(g);
        DistMatrix<T,Star,MC  > x1_Star_MC(g);
        DistMatrix<T,Star,MR  > z0_Star_MR(g);

        // Start the algorithm
        LockedPartitionUpDiagonal
        ( L, LTL, LTR,
             LBL, LBR, 0 );
        PartitionLeft( x,  xL, xR, 0 );
        while( xL.Width() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            RepartitionLeft
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );

            x1_Star_MC.AlignWith( L10 );
            z0_Star_MR.AlignWith( L10 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            blas::Trsv
            ( Lower, orientation, diagonal,
              L11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_Star_MC = x1_Star_Star;
            blas::Gemv
            ( orientation, (T)-1, 
              L10.LockedLocalMatrix(), 
              x1_Star_MC.LockedLocalMatrix(),
              (T)0, z0_Star_MR.LocalMatrix() );
            x0.SumScatterUpdate( (T)1, z0_Star_MR );
            //----------------------------------------------------------------//
            x1_Star_MC.FreeAlignments();
            z0_Star_MR.FreeAlignments();

            SlideLockedPartitionUpDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            SlidePartitionLeft
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrsvLT
( Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& x );

template void elemental::blas::internal::TrsvLT
( Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrsvLT
( Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& x );

template void elemental::blas::internal::TrsvLT
( Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

