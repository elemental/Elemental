/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_TRR2_HPP
#define ELEM_BLAS_TRR2_HPP

namespace elem {

// TODO: Generalize to both left and right diagonals
// TODO: Distributed version

// A := A + alpha X Y'
template<typename T>
inline void
Trr2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& X, const Matrix<T>& Y, Matrix<T>& A, 
  bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry entry("Trr2");
    if( X.Width() != 2 || Y.Width() != 2 )
        LogicError("X and Y must be of width 2");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
#ifndef RELEASE
    if( X.Height() != m || Y.Height() != n )
        LogicError("X and Y must conform with A");
#endif
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta0 = alpha*( conjugate ? Conj(Y.Get(j,0)) : Y.Get(j,0) );
            const T eta1 = alpha*( conjugate ? Conj(Y.Get(j,1)) : Y.Get(j,1) );
            for( Int i=j; i<m; ++i )
                A.Update( i, j, X.Get(i,0)*eta0+X.Get(i,1)*eta1 );
            if( conjugate )
                A.Set( j, j, A.GetRealPart(j,j) );
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta0 = alpha*( conjugate ? Conj(Y.Get(j,0)) : Y.Get(j,0) );
            const T eta1 = alpha*( conjugate ? Conj(Y.Get(j,1)) : Y.Get(j,1) );
            for( Int i=0; i<=j; ++i )
                A.Update( i, j, X.Get(i,0)*eta0+X.Get(i,1)*eta1 );
            if( conjugate )
                A.Set( j, j, A.GetRealPart(j,j) );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_TRR2_HPP
