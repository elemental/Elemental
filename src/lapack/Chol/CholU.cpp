/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::lapack::internal::CholU
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholU");
#endif
    lapack::internal::CholUVar3( A );
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
   X11^H[MR,* ] := (A01[MC,MR])^H A01[MC,* ]
   A11[MC,MR] := A11[MC,MR] - ((SumCol(X11^H[MR,* ]))[* ,MC])^H

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X12^H[MR,* ] := (A02[MC,MR])^H A01[MC,* ]
   A12[MC,MR] := A12[MC,MR] - ((SumCol(X12^H[MR,* ]))[MC,* ])^H

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
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<T,MC,  Star> A01_MC_Star(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,Star,VR  > A12_Star_VR(g);
    DistMatrix<T,MR,  Star> X11Herm_MR_Star(g);
    DistMatrix<T,MR,  MC  > X11Herm_MR_MC(g);
    DistMatrix<T,MC,  MR  > X11(g);
    DistMatrix<T,MR,  Star> X12Herm_MR_Star(g);
    DistMatrix<T,MR,  MC  > X12Herm_MR_MC(g);
    DistMatrix<T,MC,  MR  > X12(g);

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

        A01_MC_Star.AlignWith( A01 );
        X11Herm_MR_Star.AlignWith( A01 );
        X11Herm_MR_MC.AlignWith( A11 );
        X11.AlignWith( A11 );
        X12Herm_MR_Star.AlignWith( A02 );
        X12Herm_MR_MC.AlignWith( A12 );
        X12.AlignWith( A12 );
        X11Herm_MR_Star.ResizeTo( A11.Width(), A11.Height() );
        X12Herm_MR_Star.ResizeTo( A12.Width(), A12.Height() );
        //--------------------------------------------------------------------//
        A01_MC_Star = A01;
        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (T)1, A01, A01_MC_Star, (T)0, X11Herm_MR_Star );
        X11Herm_MR_MC.SumScatterFrom( X11Herm_MR_Star );
        blas::ConjTrans( X11Herm_MR_MC, X11 );
        blas::Axpy( (T)-1, X11, A11 );

        A11_Star_Star = A11;
        lapack::internal::LocalChol( Upper, A11_Star_Star );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (T)1, A02, A01_MC_Star, (T)0, X12Herm_MR_Star );
        X12Herm_MR_MC.SumScatterFrom( X12Herm_MR_Star );
        blas::ConjTrans( X12Herm_MR_MC, X12 );
        blas::Axpy( (T)-1, X12, A12 );

        A12_Star_VR = A12;
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (T)1, A11_Star_Star, A12_Star_VR );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_MC_Star.FreeAlignments();
        X11Herm_MR_Star.FreeAlignments();
        X11Herm_MR_MC.FreeAlignments();
        X11.FreeAlignments();
        X12Herm_MR_Star.FreeAlignments();
        X12Herm_MR_MC.FreeAlignments();
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
elemental::lapack::internal::CholUVar2Naive
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholUVar2Naive");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
    if( A.GetGrid().VCRank() == 0 )
    {
        cout << "CholUVar2Naive exists solely for academic purposes. Please "
                "use CholUVar2 in real applications." << endl;
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<T,MC,  Star> A01_MC_Star(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,Star,VR  > A12_Star_VR(g);
    DistMatrix<T,Star,MR  > X11_Star_MR(g);
    DistMatrix<T,Star,MR  > X12_Star_MR(g);

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
        lapack::internal::LocalChol( Upper, A11_Star_Star );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (T)1, A01_MC_Star, A02, (T)0, X12_Star_MR );
        A12.SumScatterUpdate( (T)-1, X12_Star_MR );

        A12_Star_VR = A12;
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (T)1, A11_Star_Star, A12_Star_VR );
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

// I do not see any algorithmic optimizations to make for the upper var3 
// Cholesky, since most memory access is stride one.
template<typename T>
void
elemental::lapack::internal::CholUVar3
( DistMatrix<T,MC,MR>& A )
{ elemental::lapack::internal::CholUVar3Naive( A ); }

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
lapack::internal::CholUVar3Naive
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::CholUVar3Naive");
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

    // Temporary matrix distributions
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,Star,VR  > A12_Star_VR(g);
    DistMatrix<T,Star,MC  > A12_Star_MC(g);
    DistMatrix<T,Star,MR  > A12_Star_MR(g);

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

        A12_Star_MC.AlignWith( A22 );
        A12_Star_MR.AlignWith( A22 );
        A12_Star_VR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        lapack::internal::LocalChol( Upper, A11_Star_Star );
        A11 = A11_Star_Star;

        A12_Star_VR = A12;
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (T)1, A11_Star_Star, A12_Star_VR );

        A12_Star_MC = A12_Star_VR;
        A12_Star_MR = A12_Star_VR;
        blas::internal::LocalTriangularRankK
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

template void elemental::lapack::internal::CholUVar2Naive
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3Naive
( DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::CholU
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2Naive
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3
( DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3Naive
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::CholU
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2Naive
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3Naive
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::CholU
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar2Naive
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::CholUVar3Naive
( DistMatrix<dcomplex,MC,MR>& A );
#endif

