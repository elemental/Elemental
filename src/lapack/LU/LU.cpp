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

template<typename T>
void
elemental::lapack::LU
( DistMatrix<T,MC,MR>& A, DistMatrix<int,VC,Star>& p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
    if( A.GetGrid() != p.GetGrid() )
        throw "A and p must be distributed over the same grid.";
    if( A.Height() != p.Height() ) 
        throw "A and p must be the same height.";
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),  AB(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),  
                               A20(grid), A21(grid), A22(grid);

    DistMatrix<int,VC,Star>
        pT(grid),  p0(grid), 
        pB(grid),  p1(grid),
                   p2(grid);

    // Temporary distributions
    DistMatrix<T,  Star,Star> A11_Star_Star(grid);
    DistMatrix<T,  Star,VR  > A12_Star_VR(grid);
    DistMatrix<T,  Star,MR  > A12_Star_MR(grid);
    DistMatrix<T,  VC,  Star> A21_VC_Star(grid);
    DistMatrix<T,  Star,MC  > A21Trans_Star_MC(grid);
    DistMatrix<int,Star,Star> p1_Star_Star(grid);

    // Pivot composition
    vector<int> image;
    vector<int> preimage;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR );
    PartitionDown
    ( p, pT,
         pB );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( pT,  p0,
         /**/ /**/
               p1,
          pB,  p2 );

        AB.View1x2( ABL, ABR );

        int pivotOffset = A01.Height();
        A12_Star_VR.AlignWith( A22 );
        A12_Star_MR.AlignWith( A22 );
        A21_VC_Star.AlignWith( A22 );
        A21Trans_Star_MC.AlignWith( A22 );
        A11_Star_Star.ResizeTo( A11.Height(), A11.Width() );
        p1_Star_Star.ResizeTo( p1.Height(), 1 );
        //--------------------------------------------------------------------//
        A21_VC_Star = A21;
        A11_Star_Star = A11;

        lapack::internal::PanelLU
        ( A11_Star_Star, 
          A21_VC_Star, p1_Star_Star, pivotOffset );
        lapack::internal::ComposePivots
        ( p1_Star_Star, image, preimage, pivotOffset );
        lapack::internal::ApplyRowPivots( AB, image, preimage, pivotOffset );

        A12_Star_VR = A12;
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, Unit, (T)1, A11_Star_Star, A12_Star_VR );

        A21Trans_Star_MC.TransposeFrom( A21_VC_Star );
        A12_Star_MR = A12_Star_VR;
        blas::internal::LocalGemm
        ( Transpose, Normal, (T)-1, A21Trans_Star_MC, A12_Star_MR, (T)1, A22 );

        A11 = A11_Star_Star;
        A12 = A12_Star_MR;
        A21.TransposeFrom( A21Trans_Star_MC );
        p1 = p1_Star_Star;
        //--------------------------------------------------------------------//
        A12_Star_VR.FreeAlignments();
        A12_Star_MR.FreeAlignments();
        A21_VC_Star.FreeAlignments();
        A21Trans_Star_MC.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDown
        ( pT,  p0,
               p1,
         /**/ /**/
          pB,  p2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::LU
( DistMatrix<float,MC,MR>& A, DistMatrix<int,VC,Star>& p );

template void
elemental::lapack::LU
( DistMatrix<double,MC,MR>& A, DistMatrix<int,VC,Star>& p );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::LU
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<int,VC,Star>& p );

template void
elemental::lapack::LU
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<int,VC,Star>& p );
#endif

