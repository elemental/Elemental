/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_LVAR2_HPP
#define EL_CHOLESKY_LVAR2_HPP

#include "./LVar3.hpp"

// TODO: Reverse variants

namespace El {
namespace cholesky {

template<typename F> 
void LVar2( Matrix<F>& A )
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

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );
        auto A20 = A( ind2, ind0 );
        auto A21 = A( ind2, ind1 );

        Herk( LOWER, NORMAL, F(-1), A10, F(1), A11 );
        cholesky::LVar3Unb( A11 );
        Gemm( NORMAL, ADJOINT, F(-1), A20, A10, F(1), A21 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11, A21 );
    }
}

template<typename F> 
void LVar2( AbstractDistMatrix<F>& APre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("Can only compute Cholesky factor of square matrices");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<F,MR,  STAR> A10Adj_MR_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,MC,  STAR> X11_MC_STAR(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );
        auto A20 = A( ind2, ind0 );
        auto A21 = A( ind2, ind1 );
 
        A10Adj_MR_STAR.AlignWith( A10 );
        A10.AdjointColAllGather( A10Adj_MR_STAR );
        X11_MC_STAR.AlignWith( A10 );
        LocalGemm( NORMAL, NORMAL, F(1), A10, A10Adj_MR_STAR, X11_MC_STAR );
        AxpyContract( F(-1), X11_MC_STAR, A11 );

        A11_STAR_STAR = A11;
        Cholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        X21_MC_STAR.AlignWith( A20 );
        LocalGemm( NORMAL, NORMAL, F(1), A20, A10Adj_MR_STAR, X21_MC_STAR );
        AxpyContract( F(-1), X21_MC_STAR, A21 );

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_LVAR2_HPP
