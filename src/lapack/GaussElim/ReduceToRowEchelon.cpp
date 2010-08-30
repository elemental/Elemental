/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::lapack::internal::ReduceToRowEchelon
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::ReduceToRowEchelon");
    if( A.GetGrid() != B.GetGrid() )
        throw logic_error( "A and B must be distributed over the same grid." );
    if( A.Height() != B.Height() )
        throw logic_error( "A and B must be the same height." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  APan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR>
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);

    // Temporary distributions
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,Star,VR  > A12_Star_VR(g);
    DistMatrix<T,Star,MR  > A12_Star_MR(g);
    DistMatrix<T,VC,  Star> A21_VC_Star(g);
    DistMatrix<T,Star,MC  > A21Trans_Star_MC(g);
    DistMatrix<T,Star,VR  > B1_Star_VR(g);
    DistMatrix<T,Star,MR  > B1_Star_MR(g);
    DistMatrix<int,Star,Star> p1_Star_Star(g);

    // In case B's columns are not aligned with A's
    const bool BAligned = ( B.ColShift() == A.ColShift() );
    DistMatrix<T,Star,MC> A21Trans_Star_MC_B(g);

    // Pivot composition
    vector<int> image;
    vector<int> preimage;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( B, BT,
         BB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );

        APan.View2x1
        ( A12,
          A22 );


        A12_Star_VR.AlignWith( A22 );
        A12_Star_MR.AlignWith( A22 );
        A21_VC_Star.AlignWith( A22 );
        A21Trans_Star_MC.AlignWith( A22 );
        B1_Star_VR.AlignWith( B1 );
        B1_Star_MR.AlignWith( B1 );
        if( ! BAligned )
            A21Trans_Star_MC_B.AlignWith( B2 );
        A11_Star_Star.ResizeTo( A11.Height(), A11.Width() );
        p1_Star_Star.ResizeTo( A11.Height(), 1 );
        //--------------------------------------------------------------------//
        A21_VC_Star = A21;
        A11_Star_Star = A11;

        lapack::internal::PanelLU
        ( A11_Star_Star,
          A21_VC_Star, p1_Star_Star, A00.Height() );
        lapack::internal::ComposePivots
        ( p1_Star_Star, image, preimage, A00.Height() );
        lapack::internal::ApplyRowPivots( APan, image, preimage, A00.Height() );
        lapack::internal::ApplyRowPivots( BB,   image, preimage, A00.Height() );

        A12_Star_VR = A12;
        B1_Star_VR = B1;
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, Unit, (T)1, A11_Star_Star, A12_Star_VR );
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, Unit, (T)1, A11_Star_Star, B1_Star_VR );

        A21Trans_Star_MC.TransposeFrom( A21_VC_Star );
        A12_Star_MR = A12_Star_VR;
        B1_Star_MR = B1_Star_VR;
        blas::internal::LocalGemm
        ( Transpose, Normal, (T)-1, A21Trans_Star_MC, A12_Star_MR, (T)1, A22 );
        if( BAligned )
        {
            blas::internal::LocalGemm
            ( Transpose, Normal, 
              (T)-1, A21Trans_Star_MC, B1_Star_MR, (T)1, B2 );
        }
        else
        {
            A21Trans_Star_MC_B = A21Trans_Star_MC;
            blas::internal::LocalGemm
            ( Transpose, Normal, 
              (T)-1, A21Trans_Star_MC_B, B1_Star_MR, (T)1, B2 );
        }

        A11 = A11_Star_Star;
        A12 = A12_Star_MR;
        B1 = B1_Star_MR;
        //--------------------------------------------------------------------//
        A12_Star_VR.FreeAlignments();
        A12_Star_MR.FreeAlignments();
        A21_VC_Star.FreeAlignments();
        A21Trans_Star_MC.FreeAlignments();
        B1_Star_VR.FreeAlignments();
        B1_Star_MR.FreeAlignments();
        if( ! BAligned )
            A21Trans_Star_MC_B.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDown
        ( BT,  B0,
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

