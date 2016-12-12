/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_REVERSE_LOWER_VARIANT3_HPP
#define EL_CHOLESKY_REVERSE_LOWER_VARIANT3_HPP

namespace El {
namespace cholesky {

template<typename F>
void ReverseLowerVariant3Unblocked( Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

        for( Int k=0; k<j; ++k )
            A(j,k) /= alpha;

        for( Int k=0; k<j; ++k )
            for( Int i=k; i<j; ++i )
                A(i,k) -= A(j,k)*Conj(A(j,i));
    }
}

template<typename F>
void ReverseLowerVariant3Blocked( Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

        cholesky::ReverseLowerVariant3Unblocked( A11 );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A11, A10 );
        Herk( LOWER, ADJOINT, Base<F>(-1), A10, Base<F>(1), A00 );
    }
}

template<typename F>
void ReverseLowerVariant3Blocked( AbstractDistMatrix<F>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )
    const Grid& grid = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(grid);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(grid);
    DistMatrix<F,STAR,MC  > A10_STAR_MC(grid);
    DistMatrix<F,STAR,MR  > A10_STAR_MR(grid);

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

#endif // ifndef EL_CHOLESKY_REVERSE_LOWER_VARIANT3_HPP
