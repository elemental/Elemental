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
TrdtrmmUVar1( Orientation orientation, Matrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TrtdrmmUVar1");
    if( U.Height() != U.Width() )
        throw std::logic_error("U must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transpose");
#endif
    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;
    Matrix<F> d1, S01;

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

        //--------------------------------------------------------------------/
        U11.GetDiagonal( d1 );
        S01 = U01;
        DiagonalSolve( LEFT, NORMAL, d1, U01, true );
        Trrk( UPPER, NORMAL, orientation, F(1), U01, S01, F(1), U00 );
        Trmm( RIGHT, UPPER, ADJOINT, UNIT, F(1), U11, U01 );
        TrdtrmmUUnblocked( orientation, U11 );
        //--------------------------------------------------------------------/

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

template<typename F>
inline void
TrdtrmmUVar1( Orientation orientation, DistMatrix<F,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TrdtrmmUVar1");
    if( U.Height() != U.Width() )
        throw std::logic_error("U must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transpose");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);
    DistMatrix<F,MD,STAR> d1(g);

    // Temporary distributions
    DistMatrix<F,MC,  STAR> S01_MC_STAR(g);
    DistMatrix<F,VC,  STAR> S01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<F,STAR,MR  > U01AdjOrTrans_STAR_MR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);

    S01_MC_STAR.AlignWith( U );
    S01_VC_STAR.AlignWith( U );
    U01_VR_STAR.AlignWith( U );
    U01AdjOrTrans_STAR_MR.AlignWith( U );

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

        //--------------------------------------------------------------------//
        U11.GetDiagonal( d1 );
        S01_MC_STAR = U01;
        S01_VC_STAR = S01_MC_STAR;
        U01_VR_STAR = S01_VC_STAR;
        if( orientation == TRANSPOSE )
        {
            DiagonalSolve( RIGHT, NORMAL, d1, U01_VR_STAR );
            U01AdjOrTrans_STAR_MR.TransposeFrom( U01_VR_STAR );
        }
        else
        {
            DiagonalSolve( RIGHT, ADJOINT, d1, U01_VR_STAR );
            U01AdjOrTrans_STAR_MR.AdjointFrom( U01_VR_STAR );
        }
        LocalTrrk( UPPER, F(1), S01_MC_STAR, U01AdjOrTrans_STAR_MR, F(1), U00 );

        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, ADJOINT, UNIT, F(1), U11_STAR_STAR, U01_VR_STAR );
        U01 = U01_VR_STAR;

        LocalTrdtrmm( orientation, UPPER, U11_STAR_STAR );
        U11 = U11_STAR_STAR;
        //--------------------------------------------------------------------//
        d1.FreeAlignments();

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

} // namespace internal
} // namespace elem
