/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CHOLESKY_UVAR3_HPP
#define ELEM_CHOLESKY_UVAR3_HPP

#include ELEM_HERK_INC
#include ELEM_TRSM_INC

namespace elem {
namespace cholesky {

template<typename F>
inline void
UVar3Unb( Matrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar3Unb");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int lda = A.LDim();
    F* ABuffer = A.Buffer();
    for( Int j=0; j<n; ++j )
    {
        Real alpha = RealPart(ABuffer[j+j*lda]);
        if( alpha <= Real(0) )
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
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::ReverseUVar3Unb");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int lda = A.LDim();
    F* ABuffer = A.Buffer();
    for( Int j=n-1; j>=0; --j )
    {
        Real alpha = RealPart(ABuffer[j+j*lda]);
        if( alpha <= Real(0) )
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
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar3");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A22 = ViewRange( A, k+nb, k+nb, n,    n    );

        cholesky::UVar3Unb( A11 );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11, A12 );
        Herk( UPPER, ADJOINT, F(-1), A12, F(1), A22 );
    }
}

template<typename F> 
inline void
ReverseUVar3( Matrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::ReverseUVar3");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto A00 = ViewRange( A, 0, 0, k,    k    );
        auto A01 = ViewRange( A, 0, k, k,    k+nb );
        auto A11 = ViewRange( A, k, k, k+nb, k+nb );

        cholesky::ReverseUVar3Unb( A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A01 );
        Herk( UPPER, NORMAL, F(-1), A01, F(1), A00 );
    }
}


template<typename F> 
inline void
UVar3( DistMatrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::UVar3");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A12 = ViewRange( A, k,    k+nb, k+nb, n    );
        auto A22 = ViewRange( A, k+nb, k+nb, n,    n    );

        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A12_STAR_MC.AlignWith( A22 );
        A12_STAR_MC = A12_STAR_VR;
        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        LocalTrrk
        ( UPPER, ADJOINT, F(-1), A12_STAR_MC, A12_STAR_MR, F(1), A22 );
        A12 = A12_STAR_MR;
    }
}

template<typename F> 
inline void
ReverseUVar3( DistMatrix<F>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::ReverseUVar3");
        if( A.Height() != A.Width() )
            LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A01_VR_STAR(g);
    DistMatrix<F,STAR,MC  > A01Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A01Adj_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto A00 = ViewRange( A, 0, 0, k,    k    );
        auto A01 = ViewRange( A, 0, k, k,    k+nb );
        auto A11 = ViewRange( A, k, k, k+nb, k+nb );

        A11_STAR_STAR = A11;
        LocalReverseCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A01_VC_STAR.AlignWith( A00 );
        A01_VC_STAR = A01;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A01_VC_STAR );

        A01_VR_STAR.AlignWith( A00 );
        A01_VR_STAR = A01_VC_STAR; 
        A01Trans_STAR_MC.AlignWith( A00 );
        A01Adj_STAR_MR.AlignWith( A00 );
        A01_VC_STAR.TransposePartialColAllGather( A01Trans_STAR_MC );
        A01_VR_STAR.AdjointPartialColAllGather( A01Adj_STAR_MR );
        LocalTrrk
        ( UPPER, TRANSPOSE, 
          F(-1), A01Trans_STAR_MC, A01Adj_STAR_MR, F(1), A00 );
        A01.TransposeRowFilterFrom( A01Trans_STAR_MC );
    }
}

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_CHOLESKY_UVAR3_HPP
