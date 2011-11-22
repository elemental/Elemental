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

template<typename T>
inline void
elemental::basic::internal::HetrmmUVar1( DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("basic::internal::HetrmmUVar1");
    if( U.Height() != U.Width() )
        throw std::logic_error("U must be square");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<T,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<T,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<T,STAR,MR  > U01Adj_STAR_MR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);

    PartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( UTL.Height() < U.Height() && UTL.Width() < U.Height() )
    {
        RepartitionDownDiagonal 
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        U01_MC_STAR.AlignWith( U00 );
        U01_VC_STAR.AlignWith( U00 );
        U01_VR_STAR.AlignWith( U00 );
        U01Adj_STAR_MR.AlignWith( U00 );
        //--------------------------------------------------------------------//
        U01_MC_STAR = U01;
        U01_VC_STAR = U01_MC_STAR;
        U01_VR_STAR = U01_VC_STAR;
        U01Adj_STAR_MR.AdjointFrom( U01_VR_STAR );
        basic::internal::LocalTrrk
        ( UPPER, (T)1, U01_MC_STAR, U01Adj_STAR_MR, (T)1, U22 );

        U11_STAR_STAR = U11;
        basic::internal::LocalTrmm
        ( RIGHT, UPPER, ADJOINT, NON_UNIT, (T)1, U11_STAR_STAR, U01_VC_STAR );
        U01 = U01_VC_STAR;

        basic::internal::LocalHetrmm( UPPER, U11_STAR_STAR );
        U11 = U11_STAR_STAR;
        //--------------------------------------------------------------------//
        U01Adj_STAR_MR.FreeAlignments();
        U01_VR_STAR.FreeAlignments();
        U01_VC_STAR.FreeAlignments();
        U01_MC_STAR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
