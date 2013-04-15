/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_GAUSSIANELIMINATION_HPP
#define LAPACK_GAUSSIANELIMINATION_HPP

#include "elemental/blas-like/level2/Trsv.hpp"
#include "elemental/lapack-like/LU.hpp"

namespace elem {

namespace internal {

template<typename F> 
inline void
ReduceToRowEchelon( Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("internal::ReduceToRowEchelon");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02,  APan,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;

    Matrix<F>
        BT,  B0,
        BB,  B1,
             B2;

    Matrix<int> p1;

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

        View2x1
        ( APan, A12,
                A22 );

        //--------------------------------------------------------------------//
        lu::Panel( APan, p1, A00.Height() );
        ComposePivots( p1, A00.Height(), image, preimage );
        ApplyRowPivots( BB, image, preimage );

        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, A12 );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11, B1 );

        Gemm( NORMAL, NORMAL, F(-1), A21, A12, F(1), A22 );
        Gemm( NORMAL, NORMAL, F(-1), A21, B1,  F(1), B2 );
        //--------------------------------------------------------------------//

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

template<typename F> 
inline void
ReduceToRowEchelon( DistMatrix<F>& A, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("internal::ReduceToRowEchelon");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  APan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<F>
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

        View2x1
        ( APan, A12,
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
        lu::Panel( A11_STAR_STAR, A21_MC_STAR, p1_STAR_STAR, A00.Height() );
        ComposePivots( p1_STAR_STAR, A00.Height(), image, preimage );
        ApplyRowPivots( APan, image, preimage );
        ApplyRowPivots( BB,   image, preimage );

        A12_STAR_VR = A12;
        B1_STAR_VR = B1;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );
        LocalTrsm( LEFT, LOWER, NORMAL, UNIT, F(1), A11_STAR_STAR, B1_STAR_VR );

        A12_STAR_MR = A12_STAR_VR;
        B1_STAR_MR = B1_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), A21_MC_STAR, A12_STAR_MR, F(1), A22 );
        if( BAligned )
        {
            LocalGemm
            ( NORMAL, NORMAL, F(-1), A21_MC_STAR, B1_STAR_MR, F(1), B2 );
        }
        else
        {
            A21_MC_STAR_B = A21_MC_STAR;
            LocalGemm
            ( NORMAL, NORMAL, F(-1), A21_MC_STAR_B, B1_STAR_MR, F(1), B2 );
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

} // namespace internal

template<typename F> 
inline void
GaussianElimination( Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("GaussianElimination");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    internal::ReduceToRowEchelon( A, B );
    if( B.Width() == 1 )
        Trsv( UPPER, NORMAL, NON_UNIT, A, B );
    else
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
GaussianElimination( DistMatrix<F>& A, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("GaussianElimination");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    internal::ReduceToRowEchelon( A, B );
    if( B.Width() == 1 )
        Trsv( UPPER, NORMAL, NON_UNIT, A, B );
    else
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_GAUSSIANELIMINATION_HPP
