/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T,typename S>
void ScaleTrapezoid( S alphaS, UpperOrLower uplo, Matrix<T>& A, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ScaleTrapezoid"))
    const Int height = A.Height();
    const Int width = A.Width();
    const Int ldim = A.LDim();
    const T alpha = T(alphaS);
    T* buffer = A.Buffer();

    if( uplo == UPPER )
    {
        EL_PARALLEL_FOR
        for( Int j=Max(0,offset-1); j<width; ++j )
        {
            const Int numRows = j-offset+1;
            for( Int i=0; i<numRows; ++i )
                buffer[i+j*ldim] *= alpha;
        }
    }
    else
    {
        EL_PARALLEL_FOR
        for( Int j=0; j<width; ++j )
        {
            const Int numZeroRows = Max(j-offset,0);
            for( Int i=numZeroRows; i<height; ++i )
                buffer[i+j*ldim] *= alpha;
        }
    }
}

template<typename T,typename S>
void
ScaleTrapezoid
( S alphaS, UpperOrLower uplo, AbstractDistMatrix<T>& A, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ScaleTrapezoid"))
    const Int height = A.Height();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const T alpha = T(alphaS);

    if( uplo == UPPER )
    {
        T* buffer = A.Buffer();
        const Int ldim = A.LDim();
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int lastRow = j-offset;
            const Int boundary = Min( lastRow+1, height );
            const Int numRows = A.LocalRowOffset(boundary);
            T* col = &buffer[jLoc*ldim];
            for( Int iLoc=0; iLoc<numRows; ++iLoc )
                col[iLoc] *= alpha;
        }
    }
    else
    {
        T* buffer = A.Buffer();
        const Int ldim = A.LDim();
        EL_PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int firstRow = Max(j-offset,0);
            const Int numZeroRows = A.LocalRowOffset(firstRow);
            T* col = &buffer[numZeroRows+jLoc*ldim];
            for( Int iLoc=0; iLoc<(localHeight-numZeroRows); ++iLoc )
                col[iLoc] *= alpha;
        }
    }
}

#define PROTO_TYPES(T,S) \
  template void ScaleTrapezoid \
  ( S alpha, UpperOrLower uplo, Matrix<T>& A, Int offset ); \
  template void ScaleTrapezoid \
  ( S alpha, UpperOrLower uplo, AbstractDistMatrix<T>& A, Int offset );

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T) 

#define PROTO_CPX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_TYPES(T,T)

PROTO_INT(Int);
PROTO_REAL(float);
PROTO_REAL(double);
PROTO_CPX(Complex<float>);
PROTO_CPX(Complex<double>);

} // namespace El
