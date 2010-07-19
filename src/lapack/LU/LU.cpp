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
elemental::lapack::LU
( DistMatrix<T,MC,MR>& A, DistMatrix<int,VC,Star>& p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
    if( A.GetGrid() != p.GetGrid() )
        throw logic_error( "A and p must be distributed over the same grid." );
    if( A.Height() != p.Height() ) 
        throw logic_error( "A and p must be the same height." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AB(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  
                         A20(g), A21(g), A22(g);

    DistMatrix<int,VC,Star>
        pT(g),  p0(g), 
        pB(g),  p1(g),
                p2(g);

    // Temporary distributions
    DistMatrix<T,  Star,Star> A11_Star_Star(g);
    DistMatrix<T,  Star,VR  > A12_Star_VR(g);
    DistMatrix<T,  Star,MR  > A12_Star_MR(g);
    DistMatrix<T,  VC,  Star> A21_VC_Star(g);
    DistMatrix<T,  Star,MC  > A21Trans_Star_MC(g);
    DistMatrix<int,Star,Star> p1_Star_Star(g);

    // Pivot composition
    vector<int> image;
    vector<int> preimage;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( p, pT,
         pB, 0 );
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

