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
elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvU");
#endif
    TrinvUVar3( diagonal, U );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvUVar3");
    if( U.Height() != U.Width() )
        throw logic_error( "Nonsquare matrices cannot be triangular." );
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions

    DistMatrix<T,VC,  Star> U01_VC_Star(g);
    DistMatrix<T,Star,Star> U11_Star_Star(g);
    DistMatrix<T,Star,VR  > U12_Star_VR(g);
    DistMatrix<T,Star,MC  > U01Trans_Star_MC(g);
    DistMatrix<T,MR,  Star> U12Trans_MR_Star(g);

    // Start the algorithm
    PartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( UBR.Height() < U.Height() )
    {
        RepartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        U01Trans_Star_MC.AlignWith( U02 );
        U12Trans_MR_Star.AlignWith( U02 );
        //--------------------------------------------------------------------//
        U11_Star_Star = U11;
        lapack::internal::LocalTrinv( Upper, diagonal, U11_Star_Star );
        U11 = U11_Star_Star;

        U01_VC_Star = U01;
        blas::internal::LocalTrmm
        ( Right, Upper, Normal, diagonal, (T)-1, U11_Star_Star, U01_VC_Star );

        // We transpose before the communication to avoid cache-thrashing
        // in the unpacking stage.
        U12Trans_MR_Star.TransposeFrom( U12 );
        U01Trans_Star_MC.TransposeFrom( U01_VC_Star );

        blas::internal::LocalGemm
        ( Transpose, Transpose, 
          (T)1, U01Trans_Star_MC, U12Trans_MR_Star, (T)1, U02 );
        U01.TransposeFrom( U01Trans_Star_MC );

        U12_Star_VR.TransposeFrom( U12Trans_MR_Star );
        blas::internal::LocalTrmm
        ( Left, Upper, Normal, diagonal, (T)1, U11_Star_Star, U12_Star_VR );
        U12 = U12_Star_VR;
        //--------------------------------------------------------------------//
        U01Trans_Star_MC.FreeAlignments();
        U12Trans_MR_Star.FreeAlignments();

        SlidePartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<double,MC,MR>& U );

template void elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& U );

template void elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& U );
#endif

