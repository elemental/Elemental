/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LU_FULL_HPP
#define EL_LU_FULL_HPP

namespace El {
namespace lu {

template<typename F>
void
Full
( Matrix<F>& A,
  Permutation& P,
  Permutation& Q )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);

    P.MakeIdentity( m );
    Q.MakeIdentity( n );
    P.ReserveSwaps( minDim );
    Q.ReserveSwaps( minDim );

    for( Int k=0; k<minDim; ++k )
    {
        const IR ind1( k ), ind2( k+1, END ), indB( k, END ), indR( k, END );

        // Find the index and value of the pivot candidate
        auto ABR = A( indB, indR );
        auto pivot = MaxAbsLoc( ABR );
        const Int iPiv = pivot.i + k;
        const Int jPiv = pivot.j + k;
        P.Swap( k, iPiv );
        Q.Swap( k, jPiv );

        RowSwap( A, k, iPiv );
        ColSwap( A, k, jPiv );

        // Now we can perform the update of the current panel
        const F alpha11 = A(k,k);
        auto a21 = A( ind2, ind1 );
        auto a12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        a21 *= alpha11Inv;
        Geru( F(-1), a21, a12, A22 );
    }
}

template<typename F>
void
Full
( AbstractDistMatrix<F>& APre,
  DistPermutation& P,
  DistPermutation& Q )
{
    EL_DEBUG_CSE
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    P.MakeIdentity( m );
    Q.MakeIdentity( n );
    P.ReserveSwaps( minDim );
    Q.ReserveSwaps( minDim );

    for( Int k=0; k<minDim; ++k )
    {
        const IR ind1( k ), ind2( k+1, END ), indB( k, END ), indR( k, END );

        // Find the index and value of the pivot candidate
        auto ABR = A( indB, indR );
        auto pivot = MaxAbsLoc( ABR );
        const Int iPiv = pivot.i + k;
        const Int jPiv = pivot.j + k;
        P.Swap( k, iPiv );
        Q.Swap( k, jPiv );

        RowSwap( A, iPiv, k );
        ColSwap( A, jPiv, k );

        // Now we can perform the update of the current panel
        const F alpha11 = A.Get(k,k);
        auto a21 = A( ind2, ind1 );
        auto a12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        a21 *= alpha11Inv;
        Geru( F(-1), a21, a12, A22 );
    }
}

} // namespace lu
} // namespace El

#endif // ifndef EL_LU_FULL_HPP
