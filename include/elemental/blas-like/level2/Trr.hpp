/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRR_HPP
#define ELEM_TRR_HPP

namespace elem {

// TODO: Generalize to both left and right diagonals

// A := A + alpha x y'
template<typename T>
inline void
Trr
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A, 
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trr");
        if( x.Width() != 1 || y.Width() != 1 )
            LogicError("x and y must be of width 1");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
        if( x.Height() != m || y.Height() != n )
            LogicError("x and y must conform with A");
    )
    const T* xCol = x.LockedBuffer();
    const T* yCol = y.LockedBuffer();
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta = alpha*( conjugate ? Conj(yCol[j]) : yCol[j] );
            T* ACol = A.Buffer(0,j);
            for( Int i=j; i<m; ++i )
                ACol[i] += xCol[i]*eta;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T eta = alpha*( conjugate ? Conj(yCol[j]) : yCol[j] );
            T* ACol = A.Buffer(0,j);
            for( Int i=0; i<=j; ++i )
                ACol[i] += xCol[i]*eta;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
}

template<typename T>
inline void
Trr
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A, 
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trr");
        if( x.Width() != 1 || y.Width() != 1 )
            LogicError("x and y must be of width 1");
    )
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    DEBUG_ONLY(
        if( x.Height() != A.Height() || y.Height() != A.Width() )
            LogicError("x and y must conform with A");
    )
    const T* xLocCol = x.LockedBuffer();
    const T* yLocCol = y.LockedBuffer();
    if( uplo == LOWER )
    {
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            const T eta = alpha*( conjugate ? Conj(yLocCol[j]) : yLocCol[j] );
            T* ALocCol = A.Buffer(0,jLoc);
            const Int mLocBefore = Length( j, colShift, colStride );
            for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                ALocCol[iLoc] += xLocCol[iLoc]*eta;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
    else
    {
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            const T eta = alpha*( conjugate ? Conj(yLocCol[j]) : yLocCol[j] );
            T* ALocCol = A.Buffer(0,jLoc);
            const Int mLocBefore = Length( j+1, colShift, colStride );
            for( Int iLoc=0; iLoc<mLocBefore; ++iLoc )
                ALocCol[iLoc] += xLocCol[iLoc]*eta;
            if( conjugate )
                A.MakeReal( j, j );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_TRR_HPP
