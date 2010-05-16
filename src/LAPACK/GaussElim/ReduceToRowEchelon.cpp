/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::blas;
using namespace elemental::lapack::internal;

template<typename T>
void
elemental::lapack::internal::ReduceToRowEchelon
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::ReduceToRowEchelon");
    if( A.GetGrid() != B.GetGrid() )
        throw "A and B must be distributed over the same grid.";
    if( A.Height() != B.Height() )
        throw "A and B must be the same height.";
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),  APan(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),
                               A20(grid), A21(grid), A22(grid);

    DistMatrix<T,MC,MR>
        BT(grid),  B0(grid),
        BB(grid),  B1(grid),
                   B2(grid);

    // Temporary distributions
    DistMatrix<T,Star,Star> A11_Star_Star(grid);
    DistMatrix<T,Star,VR  > A12_Star_VR(grid);
    DistMatrix<T,Star,MR  > A12_Star_MR(grid);
    DistMatrix<T,VC,  Star> A21_VC_Star(grid);
    DistMatrix<T,Star,MC  > A21Trans_Star_MC(grid);
    DistMatrix<T,Star,VR  > B1_Star_VR(grid);
    DistMatrix<T,Star,MR  > B1_Star_MR(grid);
    DistMatrix<int,Star,Star> p1_Star_Star(grid);

    // In case B's columns are not aligned with A's
    const bool BAligned = ( B.ColShift() == A.ColShift() );
    DistMatrix<T,Star,MC> A21Trans_Star_MC_B(grid);

    // Pivot composition
    vector<int> image;
    vector<int> preimage;

    // Start the algorithm
    PartitionDownDiagonal( A, ATL, ATR,
                              ABL, ABR );
    PartitionDown( B, BT,
                      BB );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                                /*************/ /******************/
                                      /**/       A10, /**/ A11, A12,
                                 ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown( BT,  B0,
                        /**/ /**/
                              B1,
                         BB,  B2 );

        APan.View2x1( A12,
                      A22 );

        A11_Star_Star.ResizeTo( A11.Height(), A11.Width() );
        A12_Star_VR.AlignWith( A22 );
        A12_Star_MR.AlignWith( A22 );
        A21_VC_Star.AlignWith( A22 );
        A21Trans_Star_MC.AlignWith( A22 );
        if( ! BAligned )
            A21Trans_Star_MC_B.AlignWith( B2 );
        B1_Star_VR.AlignWith( B1 );
        B1_Star_MR.AlignWith( B1 );
        p1_Star_Star.ResizeTo( A11.Height(), 1 );
        //--------------------------------------------------------------------//
        A21_VC_Star = A21;
        A11_Star_Star = A11;

        PanelLU( A11_Star_Star,
                 A21_VC_Star, p1_Star_Star, A00.Height() );
        ComposePivots( p1_Star_Star, image, preimage, A00.Height() );
        ApplyRowPivots( APan, image, preimage, A00.Height() );
        ApplyRowPivots( BB,   image, preimage, A00.Height() );

        A12_Star_VR = A12;
        B1_Star_VR = B1;
        Trsm( Left, Lower, Normal, Unit,
              (T)1, A11_Star_Star.LockedLocalMatrix(),
                    A12_Star_VR.LocalMatrix()         );
        Trsm( Left, Lower, Normal, Unit,
              (T)1, A11_Star_Star.LockedLocalMatrix(),
                    B1_Star_VR.LocalMatrix()          );

        A21Trans_Star_MC.TransposeFrom( A21_VC_Star );
        A12_Star_MR = A12_Star_VR;
        B1_Star_MR = B1_Star_VR;
        Gemm( Transpose, Normal,
              (T)-1, A21Trans_Star_MC.LockedLocalMatrix(),
                     A12_Star_MR.LockedLocalMatrix(),
              (T) 1, A22.LocalMatrix()                    );
        if( BAligned )
        {
            Gemm( Transpose, Normal,
                  (T)-1, A21Trans_Star_MC.LockedLocalMatrix(),
                         B1_Star_MR.LockedLocalMatrix(),
                  (T) 1, B2.LocalMatrix()                     );
        }
        else
        {
            A21Trans_Star_MC_B = A21Trans_Star_MC;
            Gemm( Transpose, Normal, 
                  (T)-1, A21Trans_Star_MC_B.LockedLocalMatrix(),
                         B1_Star_MR.LockedLocalMatrix(),
                  (T) 1, B2.LocalMatrix()                       );
        }

        A11 = A11_Star_Star;
        A12 = A12_Star_MR;
        B1 = B1_Star_MR;
        //--------------------------------------------------------------------//
        A12_Star_VR.FreeConstraints();
        A12_Star_MR.FreeConstraints();
        A21_VC_Star.FreeConstraints();
        A21Trans_Star_MC.FreeConstraints();
        if( ! BAligned )
            A21Trans_Star_MC_B.FreeConstraints();
        B1_Star_VR.FreeConstraints();
        B1_Star_MR.FreeConstraints();

        SlidePartitionDownDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                                         /**/       A10, A11, /**/ A12,
                                   /*************/ /******************/
                                    ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDown( BT,  B0,
                                 B1,
                           /**/ /**/
                            BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::internal::ReduceToRowEchelon
( DistMatrix<float,MC,MR>& A, DistMatrix<float,MC,MR>& B );

template void
elemental::lapack::internal::ReduceToRowEchelon
( DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::internal::ReduceToRowEchelon
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<scomplex,MC,MR>& B );

template void
elemental::lapack::internal::ReduceToRowEchelon
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<dcomplex,MC,MR>& B );
#endif

