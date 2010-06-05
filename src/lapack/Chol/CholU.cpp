/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

// The mainline Cholesky wraps the variant 2 algorithm
template<typename T>
void
elemental::lapack::internal::CholU
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholU");
#endif
    lapack::internal::CholUVar2( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

/*
   Parallelization of Variant 2 Upper Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A01^H A01
   A11 := Chol(A11)
   A12 := A12 - A01^H A02
   A12 := triu(A11)^-H A12
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A01[MC,* ] <- A01[MC,MR]
   X11[* ,MR] := (A01[MC,* ])^H A01[MC,MR]
   A11[MC,MR] := A11[MC,MR] - (SumCol(X11[* ,MR]))[MC,* ]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X12[* ,MR] := (A01[MC,* ])^H A02[MC,MR]
   A12[MC,MR] := A12[MC,MR] - (SumCol(X12[* ,MR]))[MC,* ]

   A12[* ,VR] <- A12[MC,MR]
   A12[* ,VR] := triu(A11[* ,* ])^-H A12[* ,VR]
   A12[MC,MR] <- A12[* ,VR]
   -----------------------------------------------------
*/
template<typename T>
void
elemental::lapack::internal::CholUVar2
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholUVar2");
    if( A.Height() != A.Width() )
        throw "Can only compute Cholesky factor of square matrices.";
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),
                               A20(grid), A21(grid), A22(grid);

    // Temporary distributions
    DistMatrix<T,MC,  Star> A01_MC_Star(grid);
    DistMatrix<T,Star,Star> A11_Star_Star(grid);
    DistMatrix<T,Star,VR  > A12_Star_VR(grid);
    DistMatrix<T,Star,MR  > X11_Star_MR(grid);
    DistMatrix<T,Star,MR  > X12_Star_MR(grid);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A01_MC_Star.AlignWith( A01 );
        X11_Star_MR.AlignWith( A01 );
        X12_Star_MR.AlignWith( A02 );
        X11_Star_MR.ResizeTo( A11.Height(), A11.Width() );
        X12_Star_MR.ResizeTo( A12.Height(), A12.Width() );
        //--------------------------------------------------------------------//
        A01_MC_Star = A01;
        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (T)1, A01_MC_Star, A01, (T)0, X11_Star_MR );
        A11.SumScatterUpdate( (T)-1, X11_Star_MR );

        A11_Star_Star = A11;
        lapack::Chol( Upper, A11_Star_Star.LocalMatrix() );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (T)1, A01_MC_Star, A02, (T)0, X12_Star_MR );
        A12.SumScatterUpdate( (T)-1, X12_Star_MR );

        A12_Star_VR = A12;
        blas::Trsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (T)1, A11_Star_Star.LockedLocalMatrix(), 
                A12_Star_VR.LocalMatrix() );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_MC_Star.FreeAlignments();
        X11_Star_MR.FreeAlignments();
        X12_Star_MR.FreeAlignments();

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
   Parallelization of Variant 3 Upper Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := Chol(A11) 
   A12 := triu(A11)^-H A12
   A22 := A22 - A12^H A12
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]
   
   A12[* ,VR] <- A12[MC,MR]
   A12[* ,VR] := triu(A11[* ,* ])^-H A12[* ,VR]

   A12[* ,VC] <- A12[* ,VR]
   A12[* ,MC] <- A12[* ,VC]
   A12[* ,MR] <- A12[* ,VR]
   A22[MC,MR] := A22[MC,MR] - (A12[* ,MC])^H A12[* ,MR]
   -----------------------------------------------------
*/
template<typename T>
void
lapack::internal::CholUVar3
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholUVar3");
    if( A.Height() != A.Width() )
        throw "Can only compute Cholesky factor of square matrices.";
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),
                               A20(grid), A21(grid), A22(grid);

    // Temporary matrix distributions
    DistMatrix<T,Star,Star> A11_Star_Star(grid);
    DistMatrix<T,Star,VR  > A12_Star_VR(grid);
    DistMatrix<T,Star,MC  > A12_Star_MC(grid);
    DistMatrix<T,Star,MR  > A12_Star_MR(grid);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR ); 
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A12_Star_MC.AlignWith( A22 );
        A12_Star_MR.AlignWith( A22 );
        A12_Star_VR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        lapack::Chol( Upper, A11_Star_Star.LocalMatrix() );
        A11 = A11_Star_Star;

        A12_Star_VR = A12;
        blas::Trsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (T)1, A11_Star_Star.LockedLocalMatrix(), 
                A12_Star_VR.LocalMatrix() );

        A12_Star_MC = A12_Star_VR;
        A12_Star_MR = A12_Star_VR;
        blas::internal::TriangularRankK
        ( Upper, ConjugateTranspose,
          (T)-1, A12_Star_MC, A12_Star_MR, (T)1, A22 );
        A12 = A12_Star_MR;
        //--------------------------------------------------------------------//
        A12_Star_MC.FreeAlignments();
        A12_Star_MR.FreeAlignments();
        A12_Star_VR.FreeAlignments();

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

template void elemental::lapack::internal::CholU
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholU
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::CholU
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholU
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3
( DistMatrix<dcomplex,MC,MR>& A );
#endif

