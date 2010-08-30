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
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvL");
#endif
    TrinvLVar3( diagonal, L );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvLVar3");
    if( L.Height() != L.Width() )
        throw logic_error( "Nonsquare matrices cannot be triangular." );
#endif
    const Grid& g = L.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<T,Star,MR  > L10_Star_MR(g);
    DistMatrix<T,Star,VR  > L10_Star_VR(g);
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,MC,  Star> L21_MC_Star(g);
    DistMatrix<T,VC,  Star> L21_VC_Star(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( LTL.Height() < L.Height() )
    {
        RepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        L10_Star_MR.AlignWith( L20 );
        L21_MC_Star.AlignWith( L20 );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11;
        lapack::internal::LocalTrinv( Lower, diagonal, L11_Star_Star );
        L11 = L11_Star_Star;

        L10_Star_VR = L10;
        blas::internal::LocalTrmm
        ( Left, Lower, Normal, diagonal, (T)-1, L11_Star_Star, L10_Star_VR );

        L21_MC_Star = L21;
        L10_Star_MR = L10_Star_VR;
        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, L21_MC_Star, L10_Star_MR, (T)1, L20 );
        L10 = L10_Star_MR;

        L21_VC_Star = L21_MC_Star;
        blas::internal::LocalTrmm
        ( Right, Lower, Normal, diagonal, (T)1, L11_Star_Star, L21_VC_Star );
        L21 = L21_VC_Star;
        //--------------------------------------------------------------------//
        L10_Star_MR.FreeAlignments();
        L21_MC_Star.FreeAlignments();

        SlidePartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<double,MC,MR>& L );

template void elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<double,MC,MR>& L );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& L );

template void elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& L );
#endif

