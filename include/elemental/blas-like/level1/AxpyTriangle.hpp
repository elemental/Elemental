/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_AXPYTRIANGLE_HPP
#define ELEM_AXPYTRIANGLE_HPP

namespace elem {

template<typename T>
inline void
AxpyTriangle( UpperOrLower uplo, T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("AxpyTriangle");
        if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
            X.Height() != Y.Height() )
            LogicError("Nonconformal AxpyTriangle");
    )
    if( uplo == UPPER )
    {
        for( Int j=0; j<X.Width(); ++j )
            blas::Axpy( j+1, alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
    }
    else
    {
        const Int n = X.Height();
        for( Int j=0; j<X.Width(); ++j )
            blas::Axpy( n-j, alpha, X.LockedBuffer(j,j), 1, Y.Buffer(j,j), 1 );
    }
}

template<typename T>
inline void
AxpyTriangle
( UpperOrLower uplo, Base<T> alpha, 
  const Matrix<T>& X, Matrix<T>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }

template<typename T,Dist U,Dist V>
inline void
AxpyTriangle
( UpperOrLower uplo, T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("AxpyTriangle");
        if( X.Grid() != Y.Grid() )
            LogicError
            ("X and Y must be distributed over the same grid");
        if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
            X.Height() != Y.Height() )
            LogicError("Nonconformal AxpyTriangle");
    )
    if( X.ColAlign() == Y.ColAlign() && X.RowAlign() == Y.RowAlign() )
    {
        const Int localHeight = X.LocalHeight();
        const Int localWidth = X.LocalWidth();
        const T* XBuffer = X.LockedBuffer();
        T* YBuffer = Y.Buffer();
        const Int XLDim = X.LDim();
        const Int YLDim = Y.LDim();
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = X.GlobalCol(jLoc);
                const Int localHeightAbove = X.LocalRowOffset(j+1);
                blas::Axpy
                ( localHeightAbove, alpha, 
                  &XBuffer[jLoc*XLDim], 1, &YBuffer[jLoc*YLDim], 1 );
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = X.GlobalCol(jLoc);
                const Int localHeightAbove = X.LocalRowOffset(j);
                const Int localHeightBelow = localHeight - localHeightAbove;
                blas::Axpy
                ( localHeightBelow, alpha, 
                  &XBuffer[localHeightAbove+jLoc*XLDim], 1,
                  &YBuffer[localHeightAbove+jLoc*YLDim], 1 );
            }
        }
    }
    else
    {
        DistMatrix<T,U,V> XCopy( X.Grid() );
        XCopy.AlignWith( Y );
        XCopy = X;
        AxpyTriangle( uplo, alpha, XCopy, Y );
    }
}

template<typename T,Dist U,Dist V>
inline void
AxpyTriangle
( UpperOrLower uplo, Base<T> alpha,
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }

} // namespace elem

#endif // ifndef ELEM_AXPYTRIANGLE_HPP
