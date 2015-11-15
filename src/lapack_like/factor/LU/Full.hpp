/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LU_FULL_HPP
#define EL_LU_FULL_HPP

namespace El {
namespace lu {

template<typename F>
inline void
Full
( Matrix<F>& A,
  Matrix<Int>& rowPiv,
  Matrix<Int>& colPiv )
{
    DEBUG_ONLY(CSE cse("lu::Full"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);

    rowPiv.Resize( minDim, 1 );
    colPiv.Resize( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        const IR ind1( k ), ind2( k+1, END ), indB( k, END ), indR( k, END );

        // Find the index and value of the pivot candidate
        auto ABR = A( indB, indR );
        auto pivot = MaxAbsLoc( ABR );
        const Int iPiv = pivot.i + k;
        const Int jPiv = pivot.j + k;
        rowPiv.Set( k, 0, iPiv );
        colPiv.Set( k, 0, jPiv );

        RowSwap( A, k, iPiv );
        ColSwap( A, k, jPiv );

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

template<typename F>
inline void
Full
( ElementalMatrix<F>& APre, 
  ElementalMatrix<Int>& rowPiv,
  ElementalMatrix<Int>& colPiv )
{
    DEBUG_ONLY(
      CSE cse("lu::Full");
      AssertSameGrids( APre, rowPiv, colPiv );
    )
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Int minDim = Min(m,n);

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    rowPiv.Resize( minDim, 1 );
    colPiv.Resize( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        const IR ind1( k ), ind2( k+1, END ), indB( k, END ), indR( k, END );

        // Find the index and value of the pivot candidate
        auto ABR = A( indB, indR );
        auto pivot = MaxAbsLoc( ABR );
        const Int iPiv = pivot.i + k;
        const Int jPiv = pivot.j + k;
        rowPiv.Set( k, 0, iPiv );
        colPiv.Set( k, 0, jPiv );

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
