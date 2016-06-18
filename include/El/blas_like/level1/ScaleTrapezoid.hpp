/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_SCALETRAPEZOID_HPP
#define EL_BLAS_SCALETRAPEZOID_HPP

namespace El {

template<typename T,typename S>
void ScaleTrapezoid( S alphaS, UpperOrLower uplo, Matrix<T>& A, Int offset )
{
    DEBUG_CSE
    if( alphaS == S(1) )
        return; 
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
    DEBUG_CSE
    if( alphaS == S(1) )
        return; 
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

template<typename T,typename S>
void
ScaleTrapezoid( S alpha, UpperOrLower uplo, SparseMatrix<T>& A, Int offset )
{
    DEBUG_CSE
    if( alpha == S(1) )
        return; 
    const Int numEntries = A.NumEntries();
    const Int* sBuf = A.LockedSourceBuffer();
    const Int *tBuf = A.LockedTargetBuffer();
    T* vBuf = A.ValueBuffer();
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = sBuf[k];
        const Int j = tBuf[k];
        if( (uplo==LOWER && j-i <= offset) || (uplo==UPPER && j-i >= offset) )
            vBuf[k] *= alpha;    
    }
}

template<typename T,typename S>
void
ScaleTrapezoid( S alpha, UpperOrLower uplo, DistSparseMatrix<T>& A, Int offset )
{
    DEBUG_CSE
    if( alpha == S(1) )
        return; 
    const Int numLocalEntries = A.NumLocalEntries();
    const Int* sBuf = A.LockedSourceBuffer();
    const Int *tBuf = A.LockedTargetBuffer();
    T* vBuf = A.ValueBuffer();
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = sBuf[k];
        const Int j = tBuf[k];
        if( (uplo==LOWER && j-i <= offset) || (uplo==UPPER && j-i >= offset) )
            vBuf[k] *= alpha;    
    }
}

} // namespace El

#endif // ifndef EL_BLAS_SCALETRAPEZOID_HPP
