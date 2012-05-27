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

/*
   Parallelization of Variant 3 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := Cholesky(A11) 
   A21 := A21 tril(A11)^-H
   A22 := A22 - A21 A21^H
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Cholesky(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]
   
   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   
   A21[VR,* ] <- A21[VC,* ]
   A21^T[* ,MC] <- A21[VC,* ]
   A21^H[* ,MR] <- A21[VR,* ]
   A22[MC,MR] := A22[MC,MR] - (A21^T[* ,MC])^T A21^H[* ,MR]
   A21[MC,MR] <- A21^T[* ,MC]
   -----------------------------------------------------
*/
template<typename F>
inline void
internal::CholeskyLVar3( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyLVar3");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);

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

        A21_VR_STAR.AlignWith( A22 );
        A21_VC_STAR.AlignWith( A22 );
        A21Trans_STAR_MC.AlignWith( A22 );
        A21Adj_STAR_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        internal::LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR = A21;
        internal::LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A21_VC_STAR );

        A21_VR_STAR = A21_VC_STAR;
        A21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        A21Adj_STAR_MR.AdjointFrom( A21_VR_STAR );

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        internal::LocalTrrk
        ( LOWER, TRANSPOSE, 
          (F)-1, A21Trans_STAR_MC, A21Adj_STAR_MR, (F)1, A22 );

        A21.TransposeFrom( A21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        A21_VR_STAR.FreeAlignments();
        A21_VC_STAR.FreeAlignments();
        A21Trans_STAR_MC.FreeAlignments();
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

/*
   Naive parallelization of Variant 3 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := Cholesky(A11) 
   A21 := A21 tril(A11)^-H
   A22 := A22 - A21 A21^H
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Cholesky(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]
   
   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   
   A21[MC,* ] <- A21[VC,* ]
   A21[MR,* ] <- A21[VC,* ]
   A22[MC,MR] := A22[MC,MR] - A21[MC,* ] (A21[MR,* ])^H
   A21[MC,MR] <- A21[MC,* ]
   -----------------------------------------------------
*/
template<typename F>
inline void
internal::CholeskyLVar3Naive( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyLVar3Naive");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
    if( A.Grid().VCRank() == 0 )
    {
        std::cout 
            << "CholeskyLVar3Naive exists solely for academic purposes. Please "
               "use CholeskyLVar3 in real applications." << std::endl;
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> A21_MR_STAR(g);

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

        A21_VC_STAR.AlignWith( A22 );
        A21_MC_STAR.AlignWith( A22 );
        A21_MR_STAR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        internal::LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR = A21;
        internal::LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A21_VC_STAR );

        A21_MC_STAR = A21_VC_STAR;
        A21_MR_STAR = A21_VC_STAR;

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        internal::LocalTrrk
        ( LOWER, ADJOINT, (F)-1, A21_MC_STAR, A21_MR_STAR, (F)1, A22 );

        A21 = A21_MC_STAR;
        //--------------------------------------------------------------------//
        A21_VC_STAR.FreeAlignments();
        A21_MC_STAR.FreeAlignments();
        A21_MR_STAR.FreeAlignments();

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

} // namespace elem
