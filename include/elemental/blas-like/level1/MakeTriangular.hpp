/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MAKETRIANGULAR_HPP
#define ELEM_MAKETRIANGULAR_HPP

namespace elem {

template<typename T>
inline void
MakeTriangular( UpperOrLower uplo, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeTriangular"))
    const Int height = A.Height();
    const Int width = A.Width();
    const Int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == LOWER )
    {
        ELEM_PARALLEL_FOR
        for( Int j=1; j<width; ++j )
        {
            const Int numZeroRows = Min( j, height );
            MemZero( &buffer[j*ldim], numZeroRows );
        }
    }
    else
    {
        ELEM_PARALLEL_FOR
        for( Int j=0; j<Min(width,height); ++j )
        {
            const Int firstZeroRow = j+1;
            MemZero( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
        }
    }
}

template<typename T,Dist U,Dist V>
inline void
MakeTriangular( UpperOrLower uplo, DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeTriangular"))
    const Int height = A.Height();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

    T* buffer = A.Buffer();
    const Int ldim = A.LDim();

    if( uplo == LOWER )
    {
        ELEM_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int lastZeroRow = j-1;
            if( lastZeroRow >= 0 )
            {
                const Int boundary = Min( lastZeroRow+1, height );
                const Int numZeroRows = A.LocalRowOffset(boundary);
                MemZero( &buffer[jLoc*ldim], numZeroRows );
            }
        }
    }
    else
    {
        ELEM_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int firstZeroRow = j+1;
            const Int numNonzeroRows = A.LocalRowOffset(firstZeroRow);
            if( numNonzeroRows < localHeight )
            {
                T* col = &buffer[numNonzeroRows+jLoc*ldim];
                MemZero( col, localHeight-numNonzeroRows );
            }
        }
    }
}

template<typename T,Dist U,Dist V>
inline void
MakeTriangular( UpperOrLower uplo, BlockDistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeTriangular"))
    const Int height = A.Height();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

    T* buffer = A.Buffer();
    const Int ldim = A.LDim();

    if( uplo == LOWER )
    {
        ELEM_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int lastZeroRow = j-1;
            if( lastZeroRow >= 0 )
            {
                const Int boundary = Min( lastZeroRow+1, height );
                const Int numZeroRows = A.LocalRowOffset(boundary);
                MemZero( &buffer[jLoc*ldim], numZeroRows );
            }
        }
    }
    else
    {
        ELEM_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int firstZeroRow = j+1;
            const Int numNonzeroRows = A.LocalRowOffset(firstZeroRow);
            if( numNonzeroRows < localHeight )
            {
                T* col = &buffer[numNonzeroRows+jLoc*ldim];
                MemZero( col, localHeight-numNonzeroRows );
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_MAKETRIANGULAR_HPP
