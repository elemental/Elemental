/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_CHOLESKY_UVAR2_HPP
#define LAPACK_CHOLESKY_UVAR2_HPP

namespace elem {
namespace internal {

template<typename F> 
inline void
CholeskyUVar2( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyUVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
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
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        //--------------------------------------------------------------------//
        Herk( UPPER, ADJOINT, F(-1), A01, F(1), A11 );
        CholeskyUVar3Unb( A11 );
        Gemm( ADJOINT, NORMAL, F(-1), A02, A01, F(1), A12 );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11, A12 );
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
CholeskyUVar2( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyUVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,MR,  STAR> X11Adj_MR_STAR(g);
    DistMatrix<F,MR,  MC  > X11Adj_MR_MC(g);
    DistMatrix<F,MR,  STAR> X12Adj_MR_STAR(g);
    DistMatrix<F,MR,  MC  > X12Adj_MR_MC(g);
    DistMatrix<F> X11(g);
    DistMatrix<F> X12(g);

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

        A01_MC_STAR.AlignWith( A01 );
        X11Adj_MR_STAR.AlignWith( A01 );
        X11Adj_MR_MC.AlignWith( A11 );
        X11.AlignWith( A11 );
        X12Adj_MR_STAR.AlignWith( A02 );
        X12Adj_MR_MC.AlignWith( A12 );
        X12.AlignWith( A12 );
        X11Adj_MR_STAR.ResizeTo( A11.Width(), A11.Height() );
        X12Adj_MR_STAR.ResizeTo( A12.Width(), A12.Height() );
        //--------------------------------------------------------------------//
        A01_MC_STAR = A01;
        LocalGemm
        ( ADJOINT, NORMAL, F(1), A01, A01_MC_STAR, F(0), X11Adj_MR_STAR );
        X11Adj_MR_MC.SumScatterFrom( X11Adj_MR_STAR );
        Adjoint( X11Adj_MR_MC, X11 );
        Axpy( F(-1), X11, A11 );

        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        LocalGemm
        ( ADJOINT, NORMAL, F(1), A02, A01_MC_STAR, F(0), X12Adj_MR_STAR );
        X12Adj_MR_MC.SumScatterFrom( X12Adj_MR_STAR );
        Adjoint( X12Adj_MR_MC, X12 );
        Axpy( F(-1), X12, A12 );

        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;
        //--------------------------------------------------------------------//
        A01_MC_STAR.FreeAlignments();
        X11Adj_MR_STAR.FreeAlignments();
        X11Adj_MR_MC.FreeAlignments();
        X11.FreeAlignments();
        X12Adj_MR_STAR.FreeAlignments();
        X12Adj_MR_MC.FreeAlignments();
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

} // namespace internal
} // namespace elem

#endif // ifndef LAPACK_CHOLESKY_UVAR2_HPP
