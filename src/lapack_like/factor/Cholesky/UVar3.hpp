/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_UVAR3_HPP
#define EL_CHOLESKY_UVAR3_HPP

namespace El {
namespace cholesky {

template<typename F>
void UVar3Unb( Matrix<F>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int ALDim = A.LDim();
    for( Int j=0; j<n; ++j )
    {
        Real alpha11 = RealPart(A(j,j));
        if( alpha11 <= Real(0) )
            LogicError("A was not numerically HPD");
        alpha11 = Sqrt( alpha11 );
        A(j,j) = alpha11;

        const Int a12Width = n-(j+1);
        F* a12 = A.Buffer(j  ,j+1);
        F* A22 = A.Buffer(j+1,j+1);

        blas::Scal( a12Width, Real(1)/alpha11, a12, ALDim );

        for( Int k=0; k<a12Width; ++k )
            A(j,j+1+k) = Conj(A(j,j+1+k));
        blas::Her( 'U', a12Width, -Real(1), a12, ALDim, A22, ALDim );
        for( Int k=0; k<a12Width; ++k )
            A(j,j+1+k) = Conj(A(j,j+1+k));
    }
}

template<typename F>
void ReverseUVar3Unb( Matrix<F>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    for( Int j=n-1; j>=0; --j )
    {
        Real alpha = RealPart(A(j,j));
        if( alpha <= Real(0) )
            LogicError("A was not numerically HPD");
        alpha = Sqrt( alpha );
        A(j,j) = alpha;
        
        // TODO: Switch to BLAS calls

        for( Int i=0; i<j; ++i )
            A(i,j) /= alpha;

        for( Int i=0; i<j; ++i )
            for( Int k=i; k<j; ++k )
                A(i,k) -= Conj(A(k,j))*A(i,j);
    }
}

template<typename F> 
void UVar3( Matrix<F>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
        auto A12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );

        cholesky::UVar3Unb( A11 );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11, A12 );
        Herk( UPPER, ADJOINT, Base<F>(-1), A12, Base<F>(1), A22 );
    }
}

template<typename F> 
void ReverseUVar3( Matrix<F>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
        auto A01 = A( ind0, ind1 );
        auto A11 = A( ind1, ind1 );

        cholesky::ReverseUVar3Unb( A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A01 );
        Herk( UPPER, NORMAL, Base<F>(-1), A01, Base<F>(1), A00 );
    }
}

template<typename F> 
void UVar3( AbstractDistMatrix<F>& APre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );

        A11_STAR_STAR = A11;
        Cholesky( UPPER, A11_STAR_STAR );
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
void ReverseUVar3( AbstractDistMatrix<F>& APre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

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

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A01 = A( ind0, ind1 );
        auto A11 = A( ind1, ind1 );

        A11_STAR_STAR = A11;
        ReverseCholesky( UPPER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        A01_VC_STAR.AlignWith( A00 );
        A01_VC_STAR = A01;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A01_VC_STAR );

        A01_VR_STAR.AlignWith( A00 );
        A01_VR_STAR = A01_VC_STAR; 
        A01Trans_STAR_MC.AlignWith( A00 );
        A01Adj_STAR_MR.AlignWith( A00 );
        Transpose( A01_VC_STAR, A01Trans_STAR_MC );
        Adjoint( A01_VR_STAR, A01Adj_STAR_MR );
        LocalTrrk
        ( UPPER, TRANSPOSE, 
          F(-1), A01Trans_STAR_MC, A01Adj_STAR_MR, F(1), A00 );
        Transpose( A01Trans_STAR_MC, A01 );
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_UVAR3_HPP
