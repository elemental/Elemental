/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {
namespace internal {

template<typename F>
inline void
TriangularInverseUVar3
( UnitOrNonUnit diag, DistMatrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TriangularInverseUVar3");
    if( U.Height() != U.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions

    DistMatrix<F,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > U12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > U01Trans_STAR_MC(g);
    DistMatrix<F,MR,  STAR> U12Trans_MR_STAR(g);

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

        U01Trans_STAR_MC.AlignWith( U02 );
        U12Trans_MR_STAR.AlignWith( U02 );
        //--------------------------------------------------------------------//
        U01_VC_STAR = U01;
        U11_STAR_STAR = U11;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, (F)-1, U11_STAR_STAR, U01_VC_STAR );

        // We transpose before the communication to avoid cache-thrashing
        // in the unpacking stage.
        U12Trans_MR_STAR.TransposeFrom( U12 );
        U01Trans_STAR_MC.TransposeFrom( U01_VC_STAR );

        LocalGemm
        ( TRANSPOSE, TRANSPOSE, 
          (F)1, U01Trans_STAR_MC, U12Trans_MR_STAR, (F)1, U02 );
        U01.TransposeFrom( U01Trans_STAR_MC );

        U12_STAR_VR.TransposeFrom( U12Trans_MR_STAR );
        LocalTrsm
        ( LEFT, UPPER, NORMAL, diag, (F)1, U11_STAR_STAR, U12_STAR_VR );
        LocalTriangularInverse( UPPER, diag, U11_STAR_STAR );
        U11 = U11_STAR_STAR;
        U12 = U12_STAR_VR;
        //--------------------------------------------------------------------//
        U01Trans_STAR_MC.FreeAlignments();
        U12Trans_MR_STAR.FreeAlignments();

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

} // namespace internal
} // namespace elem
