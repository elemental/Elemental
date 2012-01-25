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

// Start by forming the partially pivoted LU decomposition of A,
//     P A = L U,
// then inverting the system of equations,
//     inv(A) inv(P) = inv(U) inv(L),
// then,
//     inv(A) = inv(U) inv(L) P.

template<typename F> 
inline void
Inverse( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Inverse");
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot invert non-square matrices");
#endif
    Matrix<int> p;
    LU( A, p );
    TriangularInverse( UPPER, NON_UNIT, A );

    // Solve inv(A) L = inv(U) for inv(A)
    Matrix<F> ATL, ATR,
              ABL, ABR;
    Matrix<F> A00, A01, A02,
              A10, A11, A12,
              A20, A21, A22;
    Matrix<F> A1, A2;
    Matrix<F> L11,
              L21;
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height() < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        A1.View( A, 0, A00.Width(),             A.Height(), A01.Width() );
        A2.View( A, 0, A00.Width()+A01.Width(), A.Height(), A02.Width() );

        //--------------------------------------------------------------------//
        // Copy out L1
        L11 = A11;
        L21 = A21;

        // Zero the strictly lower triangular portion of A1
        A11.MakeTrapezoidal( LEFT, UPPER );
        A21.SetToZero();

        // Perform the lazy update of A1
        Gemm( NORMAL, NORMAL, (F)-1, A2, L21, (F)1, A1 );

        // Solve against this diagonal block of L11
        Trsm( RIGHT, LOWER, NORMAL, UNIT, (F)1, L11, A1 );
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /*******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }

    // inv(A) := inv(A) P
    ApplyInverseColumnPivots( A, p );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
Inverse( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("Inverse");
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot invert non-square matrices");
#endif
    const Grid& g = A.Grid();
    DistMatrix<int,VC,STAR> p( g );
    LU( A, p );
    TriangularInverse( UPPER, NON_UNIT, A );

    // Solve inv(A) L = inv(U) for inv(A)
    DistMatrix<F,MC,MR> ATL(g), ATR(g), 
                        ABL(g), ABR(g);
    DistMatrix<F,MC,MR> A00(g), A01(g), A02(g),
                        A10(g), A11(g), A12(g),
                        A20(g), A21(g), A22(g);
    DistMatrix<F,MC,MR> A1(g), A2(g);
    DistMatrix<F,VC,  STAR> A1_VC_STAR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > L21Trans_STAR_MR(g);
    DistMatrix<F,MC,  STAR> Z1(g);
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height() < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        A1.View( A, 0, A00.Width(),             A.Height(), A01.Width() );
        A2.View( A, 0, A00.Width()+A01.Width(), A.Height(), A02.Width() );

        L21_VR_STAR.AlignWith( A2 );
        L21Trans_STAR_MR.AlignWith( A2 );
        Z1.AlignWith( A01 );
        Z1.ResizeTo( A.Height(), A01.Width() );
        //--------------------------------------------------------------------//
        // Copy out L1
        L11_STAR_STAR = A11;
        L21_VR_STAR = A21;
        L21Trans_STAR_MR.TransposeFrom( L21_VR_STAR );

        // Zero the strictly lower triangular portion of A1
        A11.MakeTrapezoidal( LEFT, UPPER );
        A21.SetToZero();

        // Perform the lazy update of A1
        internal::LocalGemm
        ( NORMAL, TRANSPOSE, 
          (F)-1, A2, L21Trans_STAR_MR, (F)0, Z1 );
        A1.SumScatterUpdate( (F)1, Z1 );

        // Solve against this diagonal block of L11
        A1_VC_STAR = A1;
        internal::LocalTrsm
        ( RIGHT, LOWER, NORMAL, UNIT, (F)1, L11_STAR_STAR, A1_VC_STAR );
        A1 = A1_VC_STAR;
        //--------------------------------------------------------------------//
        Z1.FreeAlignments();
        L21Trans_STAR_MR.FreeAlignments();
        L21_VR_STAR.FreeAlignments();

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /*******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }

    // inv(A) := inv(A) P
    ApplyInverseColumnPivots( A, p );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
