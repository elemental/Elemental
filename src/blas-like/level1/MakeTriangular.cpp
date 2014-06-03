/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T>
void MakeTriangular( UpperOrLower uplo, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeTriangular"))
    const Int height = A.Height();
    const Int width = A.Width();
    const Int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == LOWER )
    {
        EL_PARALLEL_FOR
        for( Int j=1; j<width; ++j )
        {
            const Int numZeroRows = Min( j, height );
            MemZero( &buffer[j*ldim], numZeroRows );
        }
    }
    else
    {
        EL_PARALLEL_FOR
        for( Int j=0; j<Min(width,height); ++j )
        {
            const Int firstZeroRow = j+1;
            MemZero( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
        }
    }
}

template<typename T>
void MakeTriangular( UpperOrLower uplo, AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeTriangular"))
    const Int height = A.Height();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

    T* buffer = A.Buffer();
    const Int ldim = A.LDim();

    if( uplo == LOWER )
    {
        EL_PARALLEL_FOR
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
        EL_PARALLEL_FOR
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

template<typename T>
void MakeTriangular( UpperOrLower uplo, AbstractBlockDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeTriangular"))
    const Int height = A.Height();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();

    T* buffer = A.Buffer();
    const Int ldim = A.LDim();

    if( uplo == LOWER )
    {
        EL_PARALLEL_FOR
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
        EL_PARALLEL_FOR
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

#define PROTO(T) \
  template void MakeTriangular \
  ( UpperOrLower uplo, Matrix<T>& A ); \
  template void MakeTriangular \
  ( UpperOrLower uplo, AbstractDistMatrix<T>& A ); \
  template void MakeTriangular \
  ( UpperOrLower uplo, AbstractBlockDistMatrix<T>& A );

PROTO(Int);
PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
