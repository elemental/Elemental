/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename F>
inline void
CholeskyUVar3Unb( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyUVar3Unb");
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
            ABuffer[j+k*lda] /= alpha;

        for( int k=j+1; k<n; ++k )
            for( int i=j+1; i<=k; ++i )
                ABuffer[i+k*lda] -= Conj(ABuffer[j+i*lda])*ABuffer[j+k*lda];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
CholeskyUVar3( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyUVar3");
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
        CholeskyUVar3Unb( A11 );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11, A12 );
        Herk( UPPER, ADJOINT, F(-1), A12, F(1), A22 );
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
CholeskyUVar3( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::CholeskyUVar3");
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
        A12_STAR_MC.FreeAlignments();
        A12_STAR_MR.FreeAlignments();
        A12_STAR_VR.FreeAlignments();

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
