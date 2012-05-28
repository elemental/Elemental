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

// I do not see any algorithmic optimizations to make for the upper var3 
// Cholesky, since most memory access is stride one.
template<typename F> 
inline void
CholeskyUVar3( DistMatrix<F,MC,MR>& A )
{ CholeskyUVar3Naive( A ); }

/*
   Parallelization of Variant 3 Upper Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := Cholesky(A11) 
   A12 := triu(A11)^-H A12
   A22 := A22 - A12^H A12
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Cholesky(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]
   
   A12[* ,VR] <- A12[MC,MR]
   A12[* ,VR] := triu(A11[* ,* ])^-H A12[* ,VR]

   A12[* ,VC] <- A12[* ,VR]
   A12[* ,MC] <- A12[* ,VC]
   A12[* ,MR] <- A12[* ,VR]
   A22[MC,MR] := A22[MC,MR] - (A12[* ,MC])^H A12[* ,MR]
   -----------------------------------------------------
*/
template<typename F> 
inline void
CholeskyUVar3Naive( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyUVar3Naive");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ( "Can only compute Cholesky factor of square matrices." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrix distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 ); 
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A12_STAR_MC.AlignWith( A22 );
        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_VR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MC = A12_STAR_VR;
        A12_STAR_MR = A12_STAR_VR;
        LocalTrrk
        ( UPPER, ADJOINT, (F)-1, A12_STAR_MC, A12_STAR_MR, (F)1, A22 );
        A12 = A12_STAR_MR;
        //--------------------------------------------------------------------//
        A12_STAR_MC.FreeAlignments();
        A12_STAR_MR.FreeAlignments();
        A12_STAR_VR.FreeAlignments();

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
