/*
   Copyright (c) 2009-2011, Jack Poulson
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

template<typename F> // representation of a real or complex number
inline void
elemental::advanced::GaussianElimination
( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("advanced::GaussianElimination");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    advanced::internal::ReduceToRowEchelon( A, B );
    if( B.Width() == 1 )
        basic::Trsv( UPPER, NORMAL, NON_UNIT, A, B );
    else
        basic::Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // representation of a real or complex number
inline void
elemental::advanced::internal::ReduceToRowEchelon
( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::ReduceToRowEchelon");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  APan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<F,MC,MR>
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,STAR,VR  > B1_STAR_VR(g);
    DistMatrix<F,STAR,MR  > B1_STAR_MR(g);
    DistMatrix<int,STAR,STAR> p1_STAR_STAR(g);

    // In case B's columns are not aligned with A's
    const bool BAligned = ( B.ColShift() == A.ColShift() );
    DistMatrix<F,MC,STAR> A21_MC_STAR_B(g);

    // Pivot composition
    std::vector<int> image, preimage;

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

        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_MR.AlignWith( A22 );
        A21_MC_STAR.AlignWith( A22 );
        B1_STAR_VR.AlignWith( B1 );
        B1_STAR_MR.AlignWith( B1 );
        if( ! BAligned )
            A21_MC_STAR_B.AlignWith( B2 );
        A11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        p1_STAR_STAR.ResizeTo( A11.Height(), 1 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        A21_MC_STAR = A21;
        advanced::internal::PanelLU
        ( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR, A00.Height() );
        advanced::internal::ComposePanelPivots
        ( p1_STAR_STAR, A00.Height(), image, preimage );
        advanced::ApplyRowPivots( APan, image, preimage );
        advanced::ApplyRowPivots( BB,   image, preimage );

        A12_STAR_VR = A12;
        B1_STAR_VR = B1;
        basic::internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, (F)1, A11_STAR_STAR, A12_STAR_VR );
        basic::internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, (F)1, A11_STAR_STAR, B1_STAR_VR );

        A12_STAR_MR = A12_STAR_VR;
        B1_STAR_MR = B1_STAR_VR;
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, A21_MC_STAR, A12_STAR_MR, (F)1, A22 );
        if( BAligned )
        {
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, (F)-1, A21_MC_STAR, B1_STAR_MR, (F)1, B2 );
        }
        else
        {
            A21_MC_STAR_B = A21_MC_STAR;
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, (F)-1, A21_MC_STAR_B, B1_STAR_MR, (F)1, B2 );
        }

        A11 = A11_STAR_STAR;
        A12 = A12_STAR_MR;
        B1 = B1_STAR_MR;
        //--------------------------------------------------------------------//
        A12_STAR_VR.FreeAlignments();
        A12_STAR_MR.FreeAlignments();
        A21_MC_STAR.FreeAlignments();
        B1_STAR_VR.FreeAlignments();
        B1_STAR_MR.FreeAlignments();
        if( ! BAligned )
            A21_MC_STAR_B.FreeAlignments();

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
