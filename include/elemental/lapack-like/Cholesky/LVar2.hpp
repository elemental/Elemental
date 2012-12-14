/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

/*
   Parallelization of Variant 2 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A10 A10^H
   A11 := Cholesky(A11)
   A21 := A21 - A20 A10^H
   A21 := A21 tril(A11)^-H
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A10^H[MR,* ] <- A10[MC,MR]
   X11[MC,* ] := A10[MC,MR] A10^H[MR,* ]
   A11[MC,MR] := A11[MC,MR] - (SumRow(X11[MC,* ]))[* ,MR]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Cholesky(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X21[MC,* ] := A20[MC,MR] A10^H[MR,* ]
   A21[MC,MR] := A21[MC,MR] - (SumRow(X21[MC,* ]))[* ,MR]

   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   A21[MC,MR] <- A21[VC,* ]
   -----------------------------------------------------
*/
template<typename F> 
inline void
CholeskyLVar2( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyLVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F> 
        ATL(g), ATR(g),   A00(g), A01(g), A02(g),
        ABL(g), ABR(g),   A10(g), A11(g), A12(g),
                          A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,MR,  STAR> A10Adj_MR_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,MC,  STAR> X11_MC_STAR(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);

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

        A10Adj_MR_STAR.AlignWith( A10 );
        X11_MC_STAR.AlignWith( A10 );
        X21_MC_STAR.AlignWith( A20 );
        X11_MC_STAR.ResizeTo( A11.Height(), A11.Width() );
        X21_MC_STAR.ResizeTo( A21.Height(), A21.Width() );
        //--------------------------------------------------------------------//
        A10Adj_MR_STAR.AdjointFrom( A10 );
        LocalGemm
        ( NORMAL, NORMAL, F(1), A10, A10Adj_MR_STAR, F(0), X11_MC_STAR );
        A11.SumScatterUpdate( F(-1), X11_MC_STAR );

        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        LocalGemm
        ( NORMAL, NORMAL, F(1), A20, A10Adj_MR_STAR, F(0), X21_MC_STAR );
        A21.SumScatterUpdate( F(-1), X21_MC_STAR );

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A10Adj_MR_STAR.FreeAlignments();
        X11_MC_STAR.FreeAlignments();
        X21_MC_STAR.FreeAlignments();

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
   Naive parallelization of Variant 2 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A10 A10^H
   A11 := Cholesky(A11)
   A21 := A21 - A20 A10^H
   A21 := A21 tril(A11)^-H
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A10[* ,MR] <- A10[MC,MR]
   X11[MC,* ] := A10[MC,MR] (A10[* ,MR])^H
   A11[MC,MR] := A11[MC,MR] - (SumRow(X11[MC,* ]))[* ,MR]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Cholesky(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X21[MC,* ] := A20[MC,MR] (A10[* ,MR])^H
   A21[MC,MR] := A21[MC,MR] - (SumRow(X21[MC,* ]))[* ,MR]

   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   A21[MC,MR] <- A21[VC,* ]
   -----------------------------------------------------
*/
template<typename F>
inline void
CholeskyLVar2Naive( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyLVar2Naive");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
    if( A.Grid().VCRank() == 0 )
    {
        std::cout 
            << "CholeskyLVar2Naive exists solely for academic purposes. Please "
               "use CholeskyLVar2 in real applications." << std::endl;
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F> 
        ATL(g), ATR(g),   A00(g), A01(g), A02(g),
        ABL(g), ABR(g),   A10(g), A11(g), A12(g),
                          A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,MC,  STAR> X11_MC_STAR(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);

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

        A10_STAR_MR.AlignWith( A10 );
        X11_MC_STAR.AlignWith( A10 );
        X21_MC_STAR.AlignWith( A20 );
        X11_MC_STAR.ResizeTo( A11.Height(), A11.Width() );
        X21_MC_STAR.ResizeTo( A21.Height(), A21.Width() );
        //--------------------------------------------------------------------//
        A10_STAR_MR = A10;
        LocalGemm( NORMAL, ADJOINT, F(1), A10, A10_STAR_MR, F(0), X11_MC_STAR );
        A11.SumScatterUpdate( F(-1), X11_MC_STAR );

        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        LocalGemm( NORMAL, ADJOINT, F(1), A20, A10_STAR_MR, F(0), X21_MC_STAR );
        A21.SumScatterUpdate( F(-1), X21_MC_STAR );

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A10_STAR_MR.FreeAlignments();
        X11_MC_STAR.FreeAlignments();
        X21_MC_STAR.FreeAlignments();

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
