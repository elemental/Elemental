/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_CHOLESKY_LVAR3_HPP
#define LAPACK_CHOLESKY_LVAR3_HPP

namespace elem {
namespace internal {

template<typename F>
inline void
CholeskyLVar3Unb( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyLVar3Unb");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
#endif
    typedef typename Base<F>::type R;

    const int n = A.Height();
    const int lda = A.LDim();
    F* ABuffer = A.Buffer();
    for( int j=0; j<n; ++j )
    {
        R alpha = RealPart(ABuffer[j+j*lda]);
        if( alpha <= R(0) )
            throw std::logic_error("A was not numerically HPD");
        alpha = Sqrt( alpha );
        ABuffer[j+j*lda] = alpha;

        for( int k=j+1; k<n; ++k )
            ABuffer[k+j*lda] /= alpha;

        for( int k=j+1; k<n; ++k )
            for( int i=k; i<n; ++i )
                ABuffer[i+k*lda] -= ABuffer[i+j*lda]*Conj(ABuffer[k+j*lda]);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
CholeskyLVar3( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyLVar3");
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
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/   
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        //--------------------------------------------------------------------//
        CholeskyLVar3Unb( A11 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11, A21 );
        Herk( LOWER, NORMAL, F(-1), A21, F(1), A22 );
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
CholeskyLVar3( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyLVar3");
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

    // Temporary matrices
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);

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

        A21_VR_STAR.AlignWith( A22 );
        A21_VC_STAR.AlignWith( A22 );
        A21Trans_STAR_MC.AlignWith( A22 );
        A21Adj_STAR_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );

        A21_VR_STAR = A21_VC_STAR;
        A21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        A21Adj_STAR_MR.AdjointFrom( A21_VR_STAR );

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        LocalTrrk
        ( LOWER, TRANSPOSE, 
          F(-1), A21Trans_STAR_MC, A21Adj_STAR_MR, F(1), A22 );

        A21.TransposeFrom( A21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        A21_VR_STAR.FreeAlignments();
        A21_VC_STAR.FreeAlignments();
        A21Trans_STAR_MC.FreeAlignments();
        A21Adj_STAR_MR.FreeAlignments();

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

#endif // ifndef LAPACK_CHOLESKY_LVAR3_HPP
