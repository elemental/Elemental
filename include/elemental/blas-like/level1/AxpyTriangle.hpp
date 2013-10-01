/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_AXPYTRIANGLE_HPP
#define ELEM_BLAS_AXPYTRIANGLE_HPP

namespace elem {

template<typename T>
inline void
AxpyTriangle( UpperOrLower uplo, T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyTriangle");
    if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
        X.Height() != Y.Height() )
        LogicError("Nonconformal AxpyTriangle");
#endif
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

#ifndef SWIG
template<typename T>
inline void
AxpyTriangle
( UpperOrLower uplo, Base<T> alpha, 
  const Matrix<T>& X, Matrix<T>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }
#endif

template<typename T,Distribution U,Distribution V>
inline void
AxpyTriangle
( UpperOrLower uplo, T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    CallStackEntry entry("AxpyTriangle");
    if( X.Grid() != Y.Grid() )
        LogicError
        ("X and Y must be distributed over the same grid");
    if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
        X.Height() != Y.Height() )
        LogicError("Nonconformal AxpyTriangle");
#endif
    if( X.ColAlign() == Y.ColAlign() && X.RowAlign() == Y.RowAlign() )
    {
        const Int localHeight = X.LocalHeight();
        const Int localWidth = X.LocalWidth();
        const Int colShift = X.ColShift();
        const Int rowShift = X.RowShift();
        const Int colStride = X.ColStride();
        const Int rowStride = X.RowStride();
        const T* XBuffer = X.LockedBuffer();
        T* YBuffer = Y.Buffer();
        const Int XLDim = X.LDim();
        const Int YLDim = Y.LDim();
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = rowShift + jLoc*rowStride;        
                const Int localHeightAbove = Length( j+1, colShift, colStride );
                blas::Axpy
                ( localHeightAbove, alpha, 
                  &XBuffer[jLoc*XLDim], 1, &YBuffer[jLoc*YLDim], 1 );
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = rowShift + jLoc*rowStride;
                const Int localHeightAbove = Length( j, colShift, colStride );
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

#ifndef SWIG
template<typename T,Distribution U,Distribution V>
inline void
AxpyTriangle
( UpperOrLower uplo, Base<T> alpha,
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }
#endif

} // namespace elem

#endif // ifndef ELEM_BLAS_AXPYTRIANGLE_HPP
