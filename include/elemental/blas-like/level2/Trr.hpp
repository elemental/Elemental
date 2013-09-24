/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_TRR_HPP
#define ELEM_BLAS_TRR_HPP

namespace elem {

// TODO: Generalize to both left and right diagonals
// TODO: Distributed version

// A := A + alpha x y'
template<typename T>
inline void
Trr
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A, 
  bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry entry("Trr");
    if( x.Width() != 1 || y.Width() != 1 )
        LogicError("x and y must be of width 1");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
#ifndef RELEASE
    if( x.Height() != m || y.Height() != n )
        LogicError("x and y must conform with A");
#endif
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta = alpha*( conjugate ? Conj(y.Get(j,0)) : y.Get(j,0) );
            for( Int i=j; i<m; ++i )
                A.Update( i, j, x.Get(i,0)*eta );
            if( conjugate )
                A.Set( j, j, A.GetRealPart(j,j) );
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta = alpha*( conjugate ? Conj(y.Get(j,0)) : y.Get(j,0) );
            for( Int i=0; i<=j; ++i )
                A.Update( i, j, x.Get(i,0)*eta );
            if( conjugate )
                A.Set( j, j, A.GetRealPart(j,j) );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_TRR_HPP
