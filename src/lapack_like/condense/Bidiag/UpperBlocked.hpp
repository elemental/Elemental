/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BIDIAG_UPPER_BLOCKED_HPP
#define EL_BIDIAG_UPPER_BLOCKED_HPP

#include "./UpperUnblocked.hpp"
#include "./UpperPanel.hpp"

namespace El {
namespace bidiag {

template<typename F>
void UpperBlocked
( Matrix<F>& A,
  Matrix<F>& householderScalarsP,
  Matrix<F>& householderScalarsQ )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    EL_DEBUG_ONLY(
      if( m < n ) 
          LogicError("A must be at least as tall as it is wide");
    )
    const Int householderScalarsPHeight = Max(n-1,0);
    const Int householderScalarsQHeight = n;
    householderScalarsP.Resize( householderScalarsPHeight, 1 );
    householderScalarsQ.Resize( householderScalarsQHeight, 1 );

    Matrix<F> X, Y;

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k, k+nb ), ind2( k+nb, END ),
                         indB( k, m    ), indR( k, END );

        auto ABR = A( indB, indR );
        auto A22 = A( ind2, ind2 );

        auto householderScalarsQ1 = householderScalarsQ( ind1, ALL );

        if( A22.Width() > 0 )
        {
            auto householderScalarsP1 = householderScalarsP( ind1, ALL );
            X.Resize( m-k, nb  );
            Y.Resize( nb,  n-k );
            bidiag::UpperPanel
            ( ABR, householderScalarsP1, householderScalarsQ1, X, Y );

            auto A12 = A( ind1, ind2 );
            auto A21 = A( ind2, ind1 );
            auto X21 = X( IR(nb,END), ALL        );
            auto Y12 = Y( ALL,        IR(nb,END) );

            // Set bottom-left entry of A12 to 1
            const F epsilon = A12(nb-1,0);
            A12(nb-1,0) = F(1);

            Gemm( NORMAL, NORMAL, F(-1), A21, Y12, F(1), A22 );
            Conjugate( A12 );
            Gemm( NORMAL, NORMAL, F(-1), X21, A12, F(1), A22 );
            Conjugate( A12 );

            // Put back bottom-left entry of A12
            A12(nb-1,0) = epsilon;
        }
        else
        {
            auto householderScalarsP1 =
              householderScalarsP( IR(k,k+nb-1), ALL );
            bidiag::UpperUnblocked
            ( ABR, householderScalarsP1, householderScalarsQ1 );
        }
    }
}

template<typename F> 
void
UpperBlocked
( DistMatrix<F>& A, 
  DistMatrix<F,STAR,STAR>& householderScalarsP,
  DistMatrix<F,STAR,STAR>& householderScalarsQ )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    EL_DEBUG_ONLY(
      AssertSameGrids( A, householderScalarsP, householderScalarsQ );
      if( m < n ) 
          LogicError("A must be at least as tall as it is wide");
    )
    const Grid& grid = A.Grid();
    const Int householderScalarsPHeight = Max(n-1,0);
    const Int householderScalarsQHeight = n;
    householderScalarsP.Resize( householderScalarsPHeight, 1 );
    householderScalarsQ.Resize( householderScalarsQHeight, 1 );
    if( grid.Size() == 1 )
    {
        UpperBlocked
        ( A.Matrix(), householderScalarsP.Matrix(),
          householderScalarsQ.Matrix() );
        return;
    }

    DistMatrix<F> X(grid), Y(grid);
    DistMatrix<F,MC,STAR> X21_MC_STAR(grid);
    DistMatrix<F,MR,STAR> Y12Adj_MR_STAR(grid);

    DistMatrix<F,MC,STAR> AB1_MC_STAR(grid);
    DistMatrix<F,MR,STAR> A1RTrans_MR_STAR(grid);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k, k+nb ), ind2( k+nb, END ),
                         indB( k, END  ), indR( k, END );

        auto A22 = A( ind2, ind2 );
        auto ABR = A( indB, indR );

        auto householderScalarsQ1 = householderScalarsQ( ind1, ALL );

        if( A22.Width() > 0 )
        {
            X.AlignWith( ABR );
            Y.AlignWith( ABR );
            X.Resize( m-k, nb  );
            Y.Resize( nb,  n-k );

            AB1_MC_STAR.AlignWith( ABR );
            A1RTrans_MR_STAR.AlignWith( ABR );
            AB1_MC_STAR.Resize( m-k, nb  );
            A1RTrans_MR_STAR.Resize( n-k, nb );

            auto householderScalarsP1 = householderScalarsP( ind1, ALL );
            bidiag::UpperPanel
            ( ABR, householderScalarsP1, householderScalarsQ1, X, Y,
              AB1_MC_STAR, A1RTrans_MR_STAR );

            auto X21 = X( IR(nb,END), ALL        );
            auto Y12 = Y( ALL,        IR(nb,END) );
            X21_MC_STAR.AlignWith( A22 );
            Y12Adj_MR_STAR.AlignWith( A22 );
            X21_MC_STAR = X21;
            Adjoint( Y12, Y12Adj_MR_STAR );

            auto A21_MC_STAR      = AB1_MC_STAR(      IR(nb,END), ALL );
            auto A12Trans_MR_STAR = A1RTrans_MR_STAR( IR(nb,END), ALL );

            LocalGemm
            ( NORMAL, ADJOINT,
              F(-1), A21_MC_STAR, Y12Adj_MR_STAR, F(1), A22 );
            LocalGemm
            ( NORMAL, ADJOINT, 
              F(-1), X21_MC_STAR, A12Trans_MR_STAR, F(1), A22 );
        }
        else
        {
            auto householderScalarsP1 =
              householderScalarsP( IR(k,k+nb-1), ALL );
            bidiag::UpperUnblocked
            ( ABR, householderScalarsP1, householderScalarsQ1 );
        }
    }
}

template<typename F> 
void
UpperBlocked
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<F>& householderScalarsPPre,
  AbstractDistMatrix<F>& householderScalarsQPre )
{
    EL_DEBUG_CSE
    DistMatrixReadWriteProxy<F,F,MC,MR>
      AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR>
      householderScalarsPProx( householderScalarsPPre ),
      householderScalarsQProx( householderScalarsQPre );
    auto& A = AProx.Get();
    auto& householderScalarsP = householderScalarsPProx.Get();
    auto& householderScalarsQ = householderScalarsQProx.Get();
    UpperBlocked( A, householderScalarsP, householderScalarsQ );
}

} // namespace bidiag
} // namespace El

#endif // ifndef EL_BIDIAG_UPPER_BLOCKED_HPP
