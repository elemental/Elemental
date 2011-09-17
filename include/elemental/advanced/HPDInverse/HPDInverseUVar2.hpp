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

template<typename F> // represents a real or complex number
inline void
elemental::advanced::internal::HPDInverseUVar2( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::HPDInverseUVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A01_VR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A01Trans_STAR_MC(g);
    DistMatrix<F,MR,  STAR> A01_MR_STAR(g);
    DistMatrix<F,STAR,MR  > A01Adj_STAR_MR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);
    DistMatrix<F,STAR,MC  > A12_STAR_MC(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A01_VC_STAR.AlignWith( A00 );
        A12_STAR_VR.AlignWith( A02 );
        A01Trans_STAR_MC.AlignWith( A00 );
        A01_VR_STAR.AlignWith( A00 );
        A01Adj_STAR_MR.AlignWith( A00 );
        A12_STAR_MR.AlignWith( A02 );
        A12_STAR_MC.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        advanced::internal::LocalCholesky( UPPER, A11_STAR_STAR );

        A01_VC_STAR = A01;
        basic::internal::LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, (F)1, A11_STAR_STAR, A01_VC_STAR );

        A12_STAR_VR = A12;
        basic::internal::LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A12_STAR_VR );

        A01Trans_STAR_MC.TransposeFrom( A01_VC_STAR );
        A01_VR_STAR = A01_VC_STAR;
        A01Adj_STAR_MR.AdjointFrom( A01_VR_STAR );
        basic::internal::LocalTriangularRankK
        ( UPPER, TRANSPOSE,
          (F)-1, A01Trans_STAR_MC, A01Adj_STAR_MR, (F)1, A00 );

        A12_STAR_MR = A12_STAR_VR;
        basic::internal::LocalGemm
        ( TRANSPOSE, NORMAL, (F)-1, A01Trans_STAR_MC, A12_STAR_MR, (F)1, A02 );

        A12_STAR_MC = A12_STAR_VR;
        basic::internal::LocalTriangularRankK
        ( UPPER, ADJOINT,
          (F)-1, A12_STAR_MC, A12_STAR_MR, (F)1, A22 );

        basic::internal::LocalTrsm
        ( RIGHT, UPPER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A01_VC_STAR );

        basic::internal::LocalTrsm
        ( LEFT, UPPER, NORMAL, NON_UNIT, (F)-1, A11_STAR_STAR, A12_STAR_VR );

        advanced::internal::LocalTriangularInverse
        ( UPPER, NON_UNIT, A11_STAR_STAR );

        basic::internal::LocalHetrmm( UPPER, A11_STAR_STAR );

        A11 = A11_STAR_STAR;
        A01 = A01_VC_STAR;
        A12 = A12_STAR_VR;
        //--------------------------------------------------------------------//
        A01_VC_STAR.FreeAlignments();
        A12_STAR_VR.FreeAlignments();
        A01Trans_STAR_MC.FreeAlignments();
        A01_VR_STAR.FreeAlignments();
        A01Adj_STAR_MR.FreeAlignments();
        A12_STAR_MR.FreeAlignments();
        A12_STAR_MC.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
