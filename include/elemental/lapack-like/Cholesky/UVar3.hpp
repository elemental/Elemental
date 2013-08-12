/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_CHOLESKY_UVAR3_HPP
#define ELEM_LAPACK_CHOLESKY_UVAR3_HPP

#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

namespace elem {
namespace cholesky {

template<typename F>
inline void
UVar3Unb( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::UVar3Unb");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    typedef BASE(F) R;

    const Int n = A.Height();
    const Int lda = A.LDim();
    F* ABuffer = A.Buffer();
    for( Int j=0; j<n; ++j )
    {
        R alpha = RealPart(ABuffer[j+j*lda]);
        if( alpha <= R(0) )
            LogicError("A was not numerically HPD");
        alpha = Sqrt( alpha );
        ABuffer[j+j*lda] = alpha;
        
        for( Int k=j+1; k<n; ++k )
            ABuffer[j+k*lda] /= alpha;

        for( Int k=j+1; k<n; ++k )
            for( Int i=j+1; i<=k; ++i )
                ABuffer[i+k*lda] -= Conj(ABuffer[j+i*lda])*ABuffer[j+k*lda];
    }
}

template<typename F>
inline void
ReverseUVar3Unb( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::ReverseUVar3Unb");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    typedef BASE(F) R;

    const Int n = A.Height();
    const Int lda = A.LDim();
    F* ABuffer = A.Buffer();
    for( Int j=n-1; j>=0; --j )
    {
        R alpha = RealPart(ABuffer[j+j*lda]);
        if( alpha <= R(0) )
            LogicError("A was not numerically HPD");
        alpha = Sqrt( alpha );
        ABuffer[j+j*lda] = alpha;
        
        for( Int i=0; i<j; ++i )
            ABuffer[i+j*lda] /= alpha;

        for( Int i=0; i<j; ++i )
            for( Int k=i; k<j; ++k )
                ABuffer[i+k*lda] -= Conj(ABuffer[k+j*lda])*ABuffer[i+j*lda];
    }
}

template<typename F> 
inline void
UVar3( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::UVar3");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
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
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        //--------------------------------------------------------------------//
        cholesky::UVar3Unb( A11 );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11, A12 );
        Herk( UPPER, ADJOINT, F(-1), A12, F(1), A22 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
}

template<typename F> 
inline void
ReverseUVar3( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::ReverseUVar3");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    // Matrix views
    Matrix<F> 
        ATL, ATR,  A00, A01, A02,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;

    // Start the algorithm
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 ); 
    while( ATL.Height() > 0 )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        //--------------------------------------------------------------------//
        cholesky::ReverseUVar3Unb( A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A01 );
        Herk( UPPER, NORMAL, F(-1), A01, F(1), A00 );
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }
}


template<typename F> 
inline void
UVar3( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::UVar3");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrix distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 ); 
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A12_STAR_MC.AlignWith( A22 );
        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_VR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MC = A12_STAR_VR;
        A12_STAR_MR = A12_STAR_VR;
        LocalTrrk
        ( UPPER, ADJOINT, F(-1), A12_STAR_MC, A12_STAR_MR, F(1), A22 );
        A12 = A12_STAR_MR;
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
}

template<typename F> 
inline void
ReverseUVar3( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::ReverseUVar3");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrix distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A01_VR_STAR(g);
    DistMatrix<F,STAR,MC  > A01Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A01Adj_STAR_MR(g);

    // Start the algorithm
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 ); 
    while( ATL.Height() > 0 )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        A01_VC_STAR.AlignWith( A00 );
        A01_VR_STAR.AlignWith( A00 );
        A01Trans_STAR_MC.AlignWith( A00 );
        A01Adj_STAR_MR.AlignWith( A00 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        LocalReverseCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A01_VC_STAR = A01;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A01_VC_STAR );

        A01_VR_STAR = A01_VC_STAR; 
        A01Trans_STAR_MC.TransposeFrom( A01_VC_STAR );
        A01Adj_STAR_MR.AdjointFrom( A01_VR_STAR );
        LocalTrrk
        ( UPPER, TRANSPOSE, 
          F(-1), A01Trans_STAR_MC, A01Adj_STAR_MR, F(1), A00 );
        A01.TransposeFrom( A01Trans_STAR_MC );
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }
}

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_LAPACK_CHOLESKY_UVAR3_HPP
