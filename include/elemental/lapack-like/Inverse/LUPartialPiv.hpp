/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_INVERSE_LUPARTIALPIV_HPP
#define LAPACK_INVERSE_LUPARTIALPIV__HPP

#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/Zero.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/ApplyColumnPivots.hpp"
#include "elemental/lapack-like/LU.hpp"
#include "elemental/lapack-like/TriangularInverse.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace inverse {

// Start by forming the partially pivoted LU decomposition of A,
//     P A = L U,
// then inverting the system of equations,
//     inv(A) inv(P) = inv(U) inv(L),
// then,
//     inv(A) = inv(U) inv(L) P.

template<typename F> 
inline void
LUPartialPiv( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("inverse::LUPartialPiv");
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot invert non-square matrices");
#endif
    Matrix<int> p;
    elem::LU( A, p );
    TriangularInverse( UPPER, NON_UNIT, A );

    // Solve inv(A) L = inv(U) for inv(A)
    Matrix<F> ATL, ATR,
              ABL, ABR;
    Matrix<F> A00, A01, A02,
              A10, A11, A12,
              A20, A21, A22;
    Matrix<F> A1, A2;
    Matrix<F> L11,
              L21;
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height() < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        View( A1, A, 0, A00.Width(),             A.Height(), A01.Width() );
        View( A2, A, 0, A00.Width()+A01.Width(), A.Height(), A02.Width() );

        //--------------------------------------------------------------------//
        // Copy out L1
        L11 = A11;
        L21 = A21;

        // Zero the strictly lower triangular portion of A1
        MakeTriangular( UPPER, A11 );
        Zero( A21 );

        // Perform the lazy update of A1
        Gemm( NORMAL, NORMAL, F(-1), A2, L21, F(1), A1 );

        // Solve against this diagonal block of L11
        Trsm( RIGHT, LOWER, NORMAL, UNIT, F(1), L11, A1 );
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /*******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }

    // inv(A) := inv(A) P
    ApplyInverseColumnPivots( A, p );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
LUPartialPiv( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("inverse::LUPartialPiv");
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot invert non-square matrices");
#endif
    const Grid& g = A.Grid();
    DistMatrix<int,VC,STAR> p( g );
    elem::LU( A, p );
    TriangularInverse( UPPER, NON_UNIT, A );

    // Solve inv(A) L = inv(U) for inv(A)
    DistMatrix<F> ATL(g), ATR(g), 
                  ABL(g), ABR(g);
    DistMatrix<F> A00(g), A01(g), A02(g),
                  A10(g), A11(g), A12(g),
                  A20(g), A21(g), A22(g);
    DistMatrix<F> A1(g), A2(g);
    DistMatrix<F,VC,  STAR> A1_VC_STAR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > L21Trans_STAR_MR(g);
    DistMatrix<F,MC,  STAR> Z1(g);
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height() < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        View( A1, A, 0, A00.Width(),             A.Height(), A01.Width() );
        View( A2, A, 0, A00.Width()+A01.Width(), A.Height(), A02.Width() );

        L21_VR_STAR.AlignWith( A2 );
        L21Trans_STAR_MR.AlignWith( A2 );
        Z1.AlignWith( A01 );
        //--------------------------------------------------------------------//
        // Copy out L1
        L11_STAR_STAR = A11;
        L21_VR_STAR = A21;
        L21Trans_STAR_MR.TransposeFrom( L21_VR_STAR );

        // Zero the strictly lower triangular portion of A1
        MakeTriangular( UPPER, A11 );
        Zero( A21 );

        // Perform the lazy update of A1
        Zeros( A.Height(), A01.Width(), Z1 );
        LocalGemm( NORMAL, TRANSPOSE, F(-1), A2, L21Trans_STAR_MR, F(0), Z1 );
        A1.SumScatterUpdate( F(1), Z1 );

        // Solve against this diagonal block of L11
        A1_VC_STAR = A1;
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, UNIT, F(1), L11_STAR_STAR, A1_VC_STAR );
        A1 = A1_VC_STAR;
        //--------------------------------------------------------------------//
        Z1.FreeAlignments();
        L21Trans_STAR_MR.FreeAlignments();
        L21_VR_STAR.FreeAlignments();

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /*******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }

    // inv(A) := inv(A) P
    ApplyInverseColumnPivots( A, p );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace inverse
} // namespace elem

#endif // ifndef LAPACK_INVERSE_LUPARTIALPIV_HPP
