/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

// The mainline Cholesky wraps the variant 3 algorithm
template<typename T>
void
elemental::lapack::internal::CholL
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholL");
#endif
    lapack::internal::CholLVar3( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

/*
   Parallelization of Variant 2 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A10 A10^H
   A11 := Chol(A11)
   A21 := A21 - A20 A10^H
   A21 := A21 tril(A11)^-H
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A10^H[MR,* ] <- A10[MC,MR]
   X11[MC,* ] := A10[MC,MR] A10^H[MR,* ]
   A11[MC,MR] := A11[MC,MR] - (SumRow(X11[MC,* ]))[* ,MR]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X21[MC,* ] := A20[MC,MR] A10^H[MR,* ]
   A21[MC,MR] := A21[MC,MR] - (SumRow(X21[MC,* ]))[* ,MR]

   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   A21[MC,MR] <- A21[VC,* ]
   -----------------------------------------------------
*/
template<typename T>
void
elemental::lapack::internal::CholLVar2
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholLVar2");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),   A00(g), A01(g), A02(g),
        ABL(g), ABR(g),   A10(g), A11(g), A12(g),
                          A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<T,MR,  Star> A10Herm_MR_Star(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,VC,  Star> A21_VC_Star(g);
    DistMatrix<T,MC,  Star> X11_MC_Star(g);
    DistMatrix<T,MC,  Star> X21_MC_Star(g);

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

        A10Herm_MR_Star.AlignWith( A10 );
        X11_MC_Star.AlignWith( A10 );
        X21_MC_Star.AlignWith( A20 );
        X11_MC_Star.ResizeTo( A11.Height(), A11.Width() );
        X21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        //--------------------------------------------------------------------//
        A10Herm_MR_Star.ConjugateTransposeFrom( A10 );
        blas::internal::LocalGemm
        ( Normal, Normal, 
          (T)1, A10, A10Herm_MR_Star, (T)0, X11_MC_Star );
        A11.SumScatterUpdate( (T)-1, X11_MC_Star );

        A11_Star_Star = A11;
        lapack::internal::LocalChol( Lower, A11_Star_Star );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( Normal, Normal,
          (T)1, A20, A10Herm_MR_Star, (T)0, X21_MC_Star );
        A21.SumScatterUpdate( (T)-1, X21_MC_Star );

        A21_VC_Star = A21;
        blas::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (T)1, A11_Star_Star, A21_VC_Star );
        A21 = A21_VC_Star;
        //--------------------------------------------------------------------//
        A10Herm_MR_Star.FreeAlignments();
        X11_MC_Star.FreeAlignments();
        X21_MC_Star.FreeAlignments();

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
   A11 := Chol(A11)
   A21 := A21 - A20 A10^H
   A21 := A21 tril(A11)^-H
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A10[* ,MR] <- A10[MC,MR]
   X11[MC,* ] := A10[MC,MR] (A10[* ,MR])^H
   A11[MC,MR] := A11[MC,MR] - (SumRow(X11[MC,* ]))[* ,MR]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X21[MC,* ] := A20[MC,MR] (A10[* ,MR])^H
   A21[MC,MR] := A21[MC,MR] - (SumRow(X21[MC,* ]))[* ,MR]

   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   A21[MC,MR] <- A21[VC,* ]
   -----------------------------------------------------
*/
template<typename T>
void
elemental::lapack::internal::CholLVar2Naive
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholLVar2Naive");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
    if( A.GetGrid().VCRank() == 0 )
    {
        cout << "CholLVar2Naive exists solely for academic purposes. Please "
                "use CholLVar2 in real applications." << endl;
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),   A00(g), A01(g), A02(g),
        ABL(g), ABR(g),   A10(g), A11(g), A12(g),
                          A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<T,Star,MR  > A10_Star_MR(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,VC,  Star> A21_VC_Star(g);
    DistMatrix<T,MC,  Star> X11_MC_Star(g);
    DistMatrix<T,MC,  Star> X21_MC_Star(g);

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

        A10_Star_MR.AlignWith( A10 );
        X11_MC_Star.AlignWith( A10 );
        X21_MC_Star.AlignWith( A20 );
        X11_MC_Star.ResizeTo( A11.Height(), A11.Width() );
        X21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        //--------------------------------------------------------------------//
        A10_Star_MR = A10;
        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          (T)1, A10, A10_Star_MR, (T)0, X11_MC_Star );
        A11.SumScatterUpdate( (T)-1, X11_MC_Star );

        A11_Star_Star = A11;
        lapack::internal::LocalChol( Lower, A11_Star_Star );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          (T)1, A20, A10_Star_MR, (T)0, X21_MC_Star );
        A21.SumScatterUpdate( (T)-1, X21_MC_Star );

        A21_VC_Star = A21;
        blas::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (T)1, A11_Star_Star, A21_VC_Star );
        A21 = A21_VC_Star;
        //--------------------------------------------------------------------//
        A10_Star_MR.FreeAlignments();
        X11_MC_Star.FreeAlignments();
        X21_MC_Star.FreeAlignments();

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
   Parallelization of Variant 3 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := Chol(A11) 
   A21 := A21 tril(A11)^-H
   A22 := A22 - A21 A21^H
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Chol(A11[* ,* ])
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
template<typename T>
void
elemental::lapack::internal::CholLVar3
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholLVar3");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrices
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,VC,  Star> A21_VC_Star(g);
    DistMatrix<T,VR,  Star> A21_VR_Star(g);
    DistMatrix<T,Star,MC  > A21Trans_Star_MC(g);
    DistMatrix<T,Star,MR  > A21Herm_Star_MR(g);

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

        A21_VR_Star.AlignWith( A22 );
        A21_VC_Star.AlignWith( A22 );
        A21Trans_Star_MC.AlignWith( A22 );
        A21Herm_Star_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        lapack::internal::LocalChol( Lower, A11_Star_Star );
        A11 = A11_Star_Star;

        A21_VC_Star = A21;
        blas::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (T)1, A11_Star_Star, A21_VC_Star );

        A21_VR_Star = A21_VC_Star;
        A21Trans_Star_MC.TransposeFrom( A21_VC_Star );
        A21Herm_Star_MR.ConjugateTransposeFrom( A21_VR_Star );

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        blas::internal::LocalTriangularRankK
        ( Lower, Transpose, 
          (T)-1, A21Trans_Star_MC, A21Herm_Star_MR, (T)1, A22 );

        A21.TransposeFrom( A21Trans_Star_MC );
        //--------------------------------------------------------------------//
        A21_VR_Star.FreeAlignments();
        A21_VC_Star.FreeAlignments();
        A21Trans_Star_MC.FreeAlignments();
        A21Herm_Star_MR.FreeAlignments();

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
   A11 := Chol(A11) 
   A21 := A21 tril(A11)^-H
   A22 := A22 - A21 A21^H
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]
   
   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   
   A21[MC,* ] <- A21[VC,* ]
   A21[MR,* ] <- A21[VC,* ]
   A22[MC,MR] := A22[MC,MR] - A21[MC,* ] (A21[MR,* ])^H
   A21[MC,MR] <- A21[MC,* ]
   -----------------------------------------------------
*/
template<typename T>
void
elemental::lapack::internal::CholLVar3Naive
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholLVar3Naive");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
    if( A.GetGrid().VCRank() == 0 )
    {
        cout << "CholLVar3Naive exists solely for academic purposes. Please "
                "use CholLVar3 in real applications." << endl;
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrices
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,VC,  Star> A21_VC_Star(g);
    DistMatrix<T,MC,  Star> A21_MC_Star(g);
    DistMatrix<T,MR,  Star> A21_MR_Star(g);

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

        A21_VC_Star.AlignWith( A22 );
        A21_MC_Star.AlignWith( A22 );
        A21_MR_Star.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        lapack::internal::LocalChol( Lower, A11_Star_Star );
        A11 = A11_Star_Star;

        A21_VC_Star = A21;
        blas::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (T)1, A11_Star_Star, A21_VC_Star );

        A21_MC_Star = A21_VC_Star;
        A21_MR_Star = A21_VC_Star;

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        blas::internal::LocalTriangularRankK
        ( Lower, ConjugateTranspose, 
          (T)-1, A21_MC_Star, A21_MR_Star, (T)1, A22 );

        A21 = A21_MC_Star;
        //--------------------------------------------------------------------//
        A21_VC_Star.FreeAlignments();
        A21_MC_Star.FreeAlignments();
        A21_MR_Star.FreeAlignments();

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

template void elemental::lapack::internal::CholL
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholLVar2
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholLVar2Naive
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholLVar3
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholLVar3Naive
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholL
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholLVar2
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholLVar2Naive
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholLVar3
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholLVar3Naive
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::CholL
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholLVar2
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholLVar2Naive
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholLVar3
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholLVar3Naive
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholL
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholLVar2
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholLVar2Naive
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholLVar3
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholLVar3Naive
( DistMatrix<dcomplex,MC,MR>& A );
#endif

