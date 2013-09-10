/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_CHOLESKY_LVAR3_HPP
#define ELEM_LAPACK_CHOLESKY_LVAR3_HPP

#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

namespace elem {
namespace cholesky {

template<typename F>
inline void
LVar3Unb( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::LVar3Unb");
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
            ABuffer[k+j*lda] /= alpha;

        for( Int k=j+1; k<n; ++k )
            for( Int i=k; i<n; ++i )
                ABuffer[i+k*lda] -= ABuffer[i+j*lda]*Conj(ABuffer[k+j*lda]);
    }
}

template<typename F>
inline void
ReverseLVar3Unb( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::ReverseLVar3Unb");
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

        for( Int k=0; k<j; ++k )
            ABuffer[j+k*lda] /= alpha;

        for( Int k=0; k<j; ++k )
            for( Int i=k; i<j; ++i )
                ABuffer[i+k*lda] -= ABuffer[j+k*lda]*Conj(ABuffer[j+i*lda]);
    }
}

template<typename F>
inline void
LVar3( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::LVar3");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = std::min(bsize,n-k);
        auto A11 = View( A, k,    k,    nb,       nb       );
        auto A21 = View( A, k+nb, k,    n-(k+nb), nb       );
        auto A22 = View( A, k+nb, k+nb, n-(k+nb), n-(k+nb) );

        cholesky::LVar3Unb( A11 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11, A21 );
        Herk( LOWER, NORMAL, F(-1), A21, F(1), A22 );
    }
}

template<typename F>
inline void
ReverseLVar3( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::ReverseLVar3");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = std::min(bsize,n-k);
        auto A00 = View( A, 0,        0,        n-(k+nb), n-(k+nb) );
        auto A10 = View( A, n-(k+nb), 0,        nb,       n-(k+nb) );
        auto A11 = View( A, n-(k+nb), n-(k+nb), nb,       nb       );

        cholesky::ReverseLVar3Unb( A11 );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A11, A10 );
        Herk( LOWER, ADJOINT, F(-1), A10, F(1), A00 );
    }
} 

template<typename F>
inline void
LVar3( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::LVar3");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = std::min(bsize,n-k);
        auto A11 = View( A, k,    k,    nb,       nb       );
        auto A21 = View( A, k+nb, k,    n-(k+nb), nb       );
        auto A22 = View( A, k+nb, k+nb, n-(k+nb), n-(k+nb) );

        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );

        A21_VR_STAR.AlignWith( A22 );
        A21_VR_STAR = A21_VC_STAR;
        A21Trans_STAR_MC.AlignWith( A22 );
        A21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        A21Adj_STAR_MR.AlignWith( A22 );
        A21Adj_STAR_MR.AdjointFrom( A21_VR_STAR );

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        LocalTrrk
        ( LOWER, TRANSPOSE, 
          F(-1), A21Trans_STAR_MC, A21Adj_STAR_MR, F(1), A22 );

        A21.TransposeFrom( A21Trans_STAR_MC );
    }
} 

template<typename F>
inline void
ReverseLVar3( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::ReverseLVar3");
    if( A.Height() != A.Width() )
        LogicError("Can only compute Cholesky factor of square matrices");
#endif
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = std::min(bsize,n-k);
        auto A00 = View( A, 0,        0,        n-(k+nb), n-(k+nb) );
        auto A10 = View( A, n-(k+nb), 0,        nb,       n-(k+nb) );
        auto A11 = View( A, n-(k+nb), n-(k+nb), nb,       nb       );

        A11_STAR_STAR = A11;
        LocalReverseCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A10_STAR_VR.AlignWith( A00 );
        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A10_STAR_VR );

        A10_STAR_MC.AlignWith( A00 );
        A10_STAR_MC = A10_STAR_VR;
        A10_STAR_MR.AlignWith( A00 );
        A10_STAR_MR = A10_STAR_VR;

        // (A10[* ,MC])^H A10[* ,MR] = (A10^H A10)[MC,MR]
        LocalTrrk
        ( LOWER, ADJOINT, F(-1), A10_STAR_MC, A10_STAR_MR, F(1), A00 );

        A10 = A10_STAR_MR;
    }
} 

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_LAPACK_CHOLESKY_LVAR3_HPP
