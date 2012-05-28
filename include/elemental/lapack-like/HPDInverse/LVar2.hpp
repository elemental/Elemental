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

//
// This approach is based upon a (conjugate)-transposition of the reordered 
// Variant 2 algorithm from Fig. 9 in Bientinesi et al.'s "Families of 
// Algorithms Related to the Inversion of a Symmetric Positive Definite Matrix".
//

template<typename F> 
inline void
HPDInverseLVar2( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HPDInverseLVar2");
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
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,STAR,MC  > A10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);

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

        A10_STAR_VR.AlignWith( A00 );
        A21_VC_STAR.AlignWith( A20 );
        A10_STAR_MC.AlignWith( A00 );
        A10_STAR_MR.AlignWith( A00 );
        A21Trans_STAR_MC.AlignWith( A20 );
        A21_VR_STAR.AlignWith( A22 );
        A21Adj_STAR_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );

        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, A11_STAR_STAR, A10_STAR_VR );

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A21_VC_STAR );

        A10_STAR_MC = A10_STAR_VR;
        A10_STAR_MR = A10_STAR_VR;
        LocalTrrk
        ( LOWER, ADJOINT,
          (F)1, A10_STAR_MC, A10_STAR_MR, (F)1, A00 );

        A21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        LocalGemm
        ( TRANSPOSE, NORMAL, (F)-1, A21Trans_STAR_MC, A10_STAR_MR, (F)1, A20 );

        A21_VR_STAR = A21_VC_STAR;
        A21Adj_STAR_MR.AdjointFrom( A21_VR_STAR );
        LocalTrrk
        ( LOWER, TRANSPOSE,
          (F)-1, A21Trans_STAR_MC, A21Adj_STAR_MR, (F)1, A22 );

        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A10_STAR_VR );

        LocalTrsm
        ( RIGHT, LOWER, NORMAL, NON_UNIT, (F)-1, A11_STAR_STAR, A21_VC_STAR );

        LocalTriangularInverse( LOWER, NON_UNIT, A11_STAR_STAR );

        LocalHetrmm( LOWER, A11_STAR_STAR );

        A11 = A11_STAR_STAR;
        A10 = A10_STAR_VR;
        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A10_STAR_VR.FreeAlignments();
        A21_VC_STAR.FreeAlignments();
        A10_STAR_MC.FreeAlignments();
        A10_STAR_MR.FreeAlignments();
        A21Trans_STAR_MC.FreeAlignments();
        A21_VR_STAR.FreeAlignments();
        A21Adj_STAR_MR.FreeAlignments();

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

} // namespace internal
} // namespace elem
