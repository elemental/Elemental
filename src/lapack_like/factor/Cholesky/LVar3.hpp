/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_LVAR3_HPP
#define EL_CHOLESKY_LVAR3_HPP

namespace El {
namespace cholesky {

template<typename F>
inline void
LVar3Unb( Matrix<F>& A )
{
    DEBUG_ONLY(
      CSE cse("cholesky::LVar3Unb");
      if( A.Height() != A.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int lda = A.LDim();
    F* ABuffer = A.Buffer();
    for( Int j=0; j<n; ++j )
    {
        Real alpha11 = RealPart(ABuffer[j+j*lda]);
        if( alpha11 <= Real(0) )
            LogicError("A was not numerically HPD");
        alpha11 = Sqrt( alpha11 );
        ABuffer[j+j*lda] = alpha11;

        const Int a21Height = n-(j+1);
        F* a21 = &ABuffer[(j+1)+ j   *lda];
        F* A22 = &ABuffer[(j+1)+(j+1)*lda];

        blas::Scal( a21Height, Real(1)/alpha11, a21, 1 );
        blas::Her( 'L', a21Height, -Real(1), a21, 1, A22, lda );
    }
}

template<typename F>
inline void
ReverseLVar3Unb( Matrix<F>& A )
{
    DEBUG_ONLY(
      CSE cse("cholesky::ReverseLVar3Unb");
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

        // TODO: Switch to BLAS calls

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
    DEBUG_ONLY(
      CSE cse("cholesky::LVar3");
      if( A.Height() != A.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        
        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        cholesky::LVar3Unb( A11 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11, A21 );
        Herk( LOWER, NORMAL, Base<F>(-1), A21, Base<F>(1), A22 );
    }
}

template<typename F>
inline void
ReverseLVar3( Matrix<F>& A )
{
    DEBUG_ONLY(
      CSE cse("cholesky::ReverseLVar3");
      if( A.Height() != A.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );

        cholesky::ReverseLVar3Unb( A11 );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A11, A10 );
        Herk( LOWER, ADJOINT, Base<F>(-1), A10, Base<F>(1), A00 );
    }
} 

template<typename F>
inline void
LVar3( AbstractDistMatrix<F>& APre )
{
    DEBUG_ONLY(
      CSE cse("cholesky::LVar3");
      if( APre.Height() != APre.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        A11_STAR_STAR = A11;
        Cholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );

        A21_VR_STAR.AlignWith( A22 );
        A21_VR_STAR = A21_VC_STAR;
        A21Trans_STAR_MC.AlignWith( A22 );
        A21Adj_STAR_MR.AlignWith( A22 );
        Transpose( A21_VC_STAR, A21Trans_STAR_MC );
        Adjoint( A21_VR_STAR, A21Adj_STAR_MR );

        // (A21^T[* ,MC])^T A21^H[* ,MR] = A21[MC,* ] A21^H[* ,MR]
        //                               = (A21 A21^H)[MC,MR]
        LocalTrrk
        ( LOWER, TRANSPOSE, 
          F(-1), A21Trans_STAR_MC, A21Adj_STAR_MR, F(1), A22 );

        Transpose( A21Trans_STAR_MC, A21 );
    }
} 

template<typename F>
inline void
ReverseLVar3( AbstractDistMatrix<F>& APre )
{
    DEBUG_ONLY(
      CSE cse("cholesky::ReverseLVar3");
      if( APre.Height() != APre.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );

        A11_STAR_STAR = A11;
        ReverseCholesky( LOWER, A11_STAR_STAR );
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
} // namespace El

#endif // ifndef EL_CHOLESKY_LVAR3_HPP
