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

/*
   Parallelization of Variant 2 Upper Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A01^H A01
   A11 := Cholesky(A11)
   A12 := A12 - A01^H A02
   A12 := triu(A11)^-H A12
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A01[MC,* ] <- A01[MC,MR]
   X11^H[MR,* ] := (A01[MC,MR])^H A01[MC,* ]
   A11[MC,MR] := A11[MC,MR] - ((SumCol(X11^H[MR,* ]))[* ,MC])^H

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Cholesky(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X12^H[MR,* ] := (A02[MC,MR])^H A01[MC,* ]
   A12[MC,MR] := A12[MC,MR] - ((SumCol(X12^H[MR,* ]))[MC,* ])^H

   A12[* ,VR] <- A12[MC,MR]
   A12[* ,VR] := triu(A11[* ,* ])^-H A12[* ,VR]
   A12[MC,MR] <- A12[* ,VR]
   -----------------------------------------------------
*/
template<typename F> 
inline void
CholeskyUVar2( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyUVar2");
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

    // Temporary distributions
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,MR,  STAR> X11Adj_MR_STAR(g);
    DistMatrix<F,MR,  MC  > X11Adj_MR_MC(g);
    DistMatrix<F,MC,  MR  > X11(g);
    DistMatrix<F,MR,  STAR> X12Adj_MR_STAR(g);
    DistMatrix<F,MR,  MC  > X12Adj_MR_MC(g);
    DistMatrix<F,MC,  MR  > X12(g);

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

        A01_MC_STAR.AlignWith( A01 );
        X11Adj_MR_STAR.AlignWith( A01 );
        X11Adj_MR_MC.AlignWith( A11 );
        X11.AlignWith( A11 );
        X12Adj_MR_STAR.AlignWith( A02 );
        X12Adj_MR_MC.AlignWith( A12 );
        X12.AlignWith( A12 );
        X11Adj_MR_STAR.ResizeTo( A11.Width(), A11.Height() );
        X12Adj_MR_STAR.ResizeTo( A12.Width(), A12.Height() );
        //--------------------------------------------------------------------//
        A01_MC_STAR = A01;
        LocalGemm
        ( ADJOINT, NORMAL, 
          (F)1, A01, A01_MC_STAR, (F)0, X11Adj_MR_STAR );
        X11Adj_MR_MC.SumScatterFrom( X11Adj_MR_STAR );
        Adjoint( X11Adj_MR_MC, X11 );
        Axpy( (F)-1, X11, A11 );

        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        LocalGemm
        ( ADJOINT, NORMAL, 
          (F)1, A02, A01_MC_STAR, (F)0, X12Adj_MR_STAR );
        X12Adj_MR_MC.SumScatterFrom( X12Adj_MR_STAR );
        Adjoint( X12Adj_MR_MC, X12 );
        Axpy( (F)-1, X12, A12 );

        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;
        //--------------------------------------------------------------------//
        A01_MC_STAR.FreeAlignments();
        X11Adj_MR_STAR.FreeAlignments();
        X11Adj_MR_MC.FreeAlignments();
        X11.FreeAlignments();
        X12Adj_MR_STAR.FreeAlignments();
        X12Adj_MR_MC.FreeAlignments();
        X12.FreeAlignments();

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
   Naive parallelization of Variant 2 Upper Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A01^H A01
   A11 := Cholesky(A11)
   A12 := A12 - A01^H A02
   A12 := triu(A11)^-H A12
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A01[MC,* ] <- A01[MC,MR]
   X11[* ,MR] := (A01[MC,* ])^H A01[MC,MR]
   A11[MC,MR] := A11[MC,MR] - (SumCol(X11[* ,MR]))[MC,* ]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Cholesky(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X12[* ,MR] := (A01[MC,* ])^H A02[MC,MR]
   A12[MC,MR] := A12[MC,MR] - (SumCol(X12[* ,MR]))[MC,* ]

   A12[* ,VR] <- A12[MC,MR]
   A12[* ,VR] := triu(A11[* ,* ])^-H A12[* ,VR]
   A12[MC,MR] <- A12[* ,VR]
   -----------------------------------------------------
*/
template<typename F> 
inline void
CholeskyUVar2Naive( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyUVar2Naive");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ( "Can only compute Cholesky factor of square matrices." );
    if( A.Grid().VCRank() == 0 )
    {
        std::cout 
            << "CholeskyUVar2Naive exists solely for academic purposes. Please "
               "use CholeskyUVar2 in real applications." << std::endl;
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MR  > X11_STAR_MR(g);
    DistMatrix<F,STAR,MR  > X12_STAR_MR(g);

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

        A01_MC_STAR.AlignWith( A01 );
        X11_STAR_MR.AlignWith( A01 );
        X12_STAR_MR.AlignWith( A02 );
        X11_STAR_MR.ResizeTo( A11.Height(), A11.Width() );
        X12_STAR_MR.ResizeTo( A12.Height(), A12.Width() );
        //--------------------------------------------------------------------//
        A01_MC_STAR = A01;
        LocalGemm
        ( ADJOINT, NORMAL, 
          (F)1, A01_MC_STAR, A01, (F)0, X11_STAR_MR );
        A11.SumScatterUpdate( (F)-1, X11_STAR_MR );

        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        LocalGemm
        ( ADJOINT, NORMAL, (F)1, A01_MC_STAR, A02, (F)0, X12_STAR_MR );
        A12.SumScatterUpdate( (F)-1, X12_STAR_MR );

        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, A11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;
        //--------------------------------------------------------------------//
        A01_MC_STAR.FreeAlignments();
        X11_STAR_MR.FreeAlignments();
        X12_STAR_MR.FreeAlignments();

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
