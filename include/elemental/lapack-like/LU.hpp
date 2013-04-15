/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_LU_HPP
#define LAPACK_LU_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/ApplyRowPivots.hpp"

#include "elemental/lapack-like/LU/Local.hpp"
#include "elemental/lapack-like/LU/Panel.hpp"

#include "elemental/lapack-like/LU/SolveAfter.hpp"

namespace elem {

template<typename F>
inline void
LocalLU( DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("LocalLU");
#endif
    LU( A.Matrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

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
        lu::Unb( A11 );
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
        LocalLU( A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_MC_STAR = A21;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A21_MC_STAR );
        A21 = A21_MC_STAR;

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR = A12_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );
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
        lu::Panel( ABRL, p1, pivotOffset );
        ComposePivots( p1, pivotOffset, image, preimage );
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

        View1x2( AB, ABL, ABR );

        const int pivotOffset = A01.Height();
        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_MR.AlignWith( A22 );
        A21_MC_STAR.AlignWith( A22 );
        A11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        p1_STAR_STAR.ResizeTo( p1.Height(), 1 );
        //--------------------------------------------------------------------//
        A21_MC_STAR = A21;
        A11_STAR_STAR = A11;
        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR, pivotOffset );
        ComposePivots( p1_STAR_STAR, pivotOffset, image, preimage );
        ApplyRowPivots( AB, image, preimage );

        // Perhaps we should give up perfectly distributing this operation since
        // it's total contribution is only O(n^2)
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MR = A12_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );

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

#endif // ifndef LAPACK_LU_HPP
