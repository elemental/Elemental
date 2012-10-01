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

#include "./LU/Local.hpp"
#include "./LU/Panel.hpp"

namespace elem {

// Performs LU factorization without pivoting

template<typename F> 
inline void
LU( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LU");
#endif
    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02, 
        ABL, ABR,  A10, A11, A12,  
                   A20, A21, A22;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        //--------------------------------------------------------------------//
        internal::UnblockedLU( A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A21 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
        //--------------------------------------------------------------------//

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

template<typename F> 
inline void
LU( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LU");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_MR.AlignWith( A22 );
        A21_MC_STAR.AlignWith( A22 );
        A11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        internal::LocalLU( A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_MC_STAR = A21;
        internal::LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A21_MC_STAR );
        A21 = A21_MC_STAR;

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR = A12;
        internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR = A12_STAR_VR;
        internal::LocalGemm
        ( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );
        A12 = A12_STAR_MR;
        //--------------------------------------------------------------------//
        A12_STAR_VR.FreeAlignments();
        A12_STAR_MR.FreeAlignments();
        A21_MC_STAR.FreeAlignments();

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

// Performs LU factorization with partial pivoting

template<typename F> 
inline void
LU( Matrix<F>& A, Matrix<int>& p )
{
#ifndef RELEASE
    PushCallStack("LU");
    if( p.Viewing() && 
        (std::min(A.Height(),A.Width()) != p.Height() || p.Width() != 1) ) 
        throw std::logic_error
        ("p must be a vector of the same height as the min dimension of A.");
#endif
    if( !p.Viewing() )
        p.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02,  ABRL, ABRR,
        ABL, ABR,  A10, A11, A12,  
                   A20, A21, A22;

    Matrix<int>
        pT,  p0, 
        pB,  p1,
             p2;

    // Pivot composition
    std::vector<int> image, preimage;

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

        PartitionRight( ABR, ABRL, ABRR, A11.Width() );

        const int pivotOffset = A01.Height();
        //--------------------------------------------------------------------//
        internal::PanelLU( ABRL, p1, pivotOffset );
        internal::ComposePanelPivots( p1, pivotOffset, image, preimage );
        ApplyRowPivots( ABL, image, preimage );
        ApplyRowPivots( ABRR, image, preimage );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
        //--------------------------------------------------------------------//

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

template<typename F> 
inline void
LU( DistMatrix<F>& A, DistMatrix<int,VC,STAR>& p )
{
#ifndef RELEASE
    PushCallStack("LU");
    if( A.Grid() != p.Grid() )
        throw std::logic_error("{A,p} must be distributed over the same grid");
    if( p.Viewing() && 
        (std::min(A.Height(),A.Width()) != p.Height() || p.Width() != 1) ) 
        throw std::logic_error
        ("p must be a vector of the same height as the min dimension of A.");
#endif
    const Grid& g = A.Grid();
    if( !p.Viewing() )
        p.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AB(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  
                         A20(g), A21(g), A22(g);

    DistMatrix<int,VC,STAR>
        pT(g),  p0(g), 
        pB(g),  p1(g),
                p2(g);

    // Temporary distributions
    DistMatrix<F,  STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,  MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,  STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,  STAR,MR  > A12_STAR_MR(g);
    DistMatrix<int,STAR,STAR> p1_STAR_STAR(g);

    // Pivot composition
    std::vector<int> image, preimage;

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

        const int pivotOffset = A01.Height();
        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_MR.AlignWith( A22 );
        A21_MC_STAR.AlignWith( A22 );
        A11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        p1_STAR_STAR.ResizeTo( p1.Height(), 1 );
        //--------------------------------------------------------------------//
        A21_MC_STAR = A21;
        A11_STAR_STAR = A11;
        internal::PanelLU
        ( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR, pivotOffset );
        internal::ComposePanelPivots
        ( p1_STAR_STAR, pivotOffset, image, preimage );
        ApplyRowPivots( AB, image, preimage );

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR = A12;
        internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR = A12_STAR_VR;
        internal::LocalGemm
        ( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );

        A11 = A11_STAR_STAR;
        A12 = A12_STAR_MR;
        A21 = A21_MC_STAR;
        p1 = p1_STAR_STAR;
        //--------------------------------------------------------------------//
        A12_STAR_VR.FreeAlignments();
        A12_STAR_MR.FreeAlignments();
        A21_MC_STAR.FreeAlignments();

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

} // namespace elem
