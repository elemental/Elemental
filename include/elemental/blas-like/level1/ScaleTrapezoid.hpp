/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_SCALETRAPEZOID_HPP
#define ELEM_BLAS_SCALETRAPEZOID_HPP

namespace elem {

template<typename T>
inline void
ScaleTrapezoid
( T alpha, UpperOrLower uplo, Matrix<T>& A, 
  Int offset=0, LeftOrRight side=LEFT )
{
    DEBUG_ONLY(CallStackEntry cse("ScaleTrapezoid"))
    const Int height = A.Height();
    const Int width = A.Width();
    const Int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == UPPER )
    {
        if( side == LEFT )
        {
            PARALLEL_FOR
            for( Int j=Max(0,offset-1); j<width; ++j )
            {
                const Int numRows = j-offset+1;
                for( Int i=0; i<numRows; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
        else
        {
            PARALLEL_FOR
            for( Int j=Max(0,offset-height+width-1); j<width; ++j )
            {
                const Int numRows = j-offset+height-width+1;
                for( Int i=0; i<numRows; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
    }
    else
    {
        if( side == LEFT )
        {
            PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                const Int numZeroRows = Max(j-offset,0);
                for( Int i=numZeroRows; i<height; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
        else
        {
            PARALLEL_FOR
            for( Int j=0; j<width; ++j )
            {
                const Int numZeroRows = Max(j-offset+height-width,0);
                for( Int i=numZeroRows; i<height; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
    }
}

template<typename T,Distribution U,Distribution V>
inline void
ScaleTrapezoid
( T alpha, UpperOrLower uplo, DistMatrix<T,U,V>& A, 
  Int offset=0, LeftOrRight side=LEFT )
{
    DEBUG_ONLY(CallStackEntry cse("ScaleTrapezoid"))
    const Int height = A.Height();
    const Int width = A.Width();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();

    if( uplo == UPPER )
    {
        T* buffer = A.Buffer();
        const Int ldim = A.LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            Int j = rowShift + jLoc*rowStride;
            Int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            Int boundary = Min( lastRow+1, height );
            Int numRows = Length_( boundary, colShift, colStride );
            T* col = &buffer[jLoc*ldim];
            for( Int iLoc=0; iLoc<numRows; ++iLoc )
                col[iLoc] *= alpha;
        }
    }
    else
    {
        T* buffer = A.Buffer();
        const Int ldim = A.LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            Int j = rowShift + jLoc*rowStride;
            Int firstRow =
                ( side==LEFT ? Max(j-offset,0)
                             : Max(j-offset+height-width,0) );
            Int numZeroRows = Length_( firstRow, colShift, colStride );
            T* col = &buffer[numZeroRows+jLoc*ldim];
            for( Int iLoc=0; iLoc<(localHeight-numZeroRows); ++iLoc )
                col[iLoc] *= alpha;
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_SCALETRAPEZOID_HPP
