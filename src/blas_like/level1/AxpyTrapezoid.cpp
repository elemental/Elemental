/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alphaS, const Matrix<T>& X, Matrix<T>& Y, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("AxpyTrapezoid");
        if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
            X.Height() != Y.Height() )
            LogicError("Nonconformal AxpyTrapezoid");
    )
    const T alpha = T(alphaS);
    const Int m = X.Height();
    const Int n = X.Width();
    if( uplo == UPPER )
    {
        for( Int j=0; j<n; ++j )
        {
            const Int iSize = Max(Min(j+1-offset,m),0);
            blas::Axpy
            ( iSize, alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const Int i = Max(Min(j-offset,m),0);
            blas::Axpy( m-i, alpha, X.LockedBuffer(i,j), 1, Y.Buffer(i,j), 1 );
        }
    }
}

template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alphaS, 
  const SparseMatrix<T>& X, SparseMatrix<T>& Y, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyTrapezoid"))
    if( X.Height() != Y.Height() || X.Width() != Y.Width() )
        LogicError("X and Y must have the same dimensions");
    const T alpha = T(alphaS);
    const Int numEntries = X.NumEntries();
    Y.Reserve( Y.NumEntries()+numEntries );
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = X.Row(k);
        const Int j = X.Col(k);
        if( (uplo==UPPER && j-i >= offset) || (uplo==LOWER && j-i <= offset) )
            Y.QueueUpdate( i, j, alpha*X.Value(k) );
    }
    Y.MakeConsistent();
}

template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alphaS, 
  const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("AxpyTrapezoid");
        AssertSameGrids( X, Y );
        if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
            X.Height() != Y.Height() )
            LogicError("Nonconformal AxpyTrapezoid");
    )
    const T alpha = T(alphaS);

    const DistData XDistData = X.DistData();
    const DistData YDistData = Y.DistData();

    if( XDistData == YDistData )
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
                const Int localHeightAbove = X.LocalRowOffset(j+1-offset);
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
                const Int localHeightAbove = X.LocalRowOffset(j-offset);
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
        std::unique_ptr<AbstractDistMatrix<T>> 
          XCopy( Y.Construct(Y.Grid(),Y.Root()) );
        XCopy->AlignWith( YDistData );
        Copy( X, *XCopy );
        AxpyTrapezoid( uplo, alpha, *XCopy, Y, offset );
    }
}

template<typename T,typename S>
void AxpyTrapezoid
( UpperOrLower uplo, S alphaS, 
  const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyTrapezoid"))
    if( X.Height() != Y.Height() || X.Width() != Y.Width() )
        LogicError("X and Y must have the same dimensions");
    if( X.Comm() != Y.Comm() )
        LogicError("X and Y must have the same communicator");
    const T alpha = T(alphaS);
    const Int numLocalEntries = X.NumLocalEntries();
    const Int firstLocalRow = X.FirstLocalRow();
    Y.Reserve( Y.NumLocalEntries()+numLocalEntries );
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = X.Row(k);
        const Int j = X.Col(k);
        if( (uplo==UPPER && j-i >= offset) || (uplo==LOWER && j-i <= offset) )
            Y.QueueLocalUpdate( i-firstLocalRow, j, alpha*X.Value(k) );
    }
    Y.MakeConsistent();
}

#define PROTO_TYPES(T,S) \
  template void AxpyTrapezoid \
  ( UpperOrLower uplo, S alpha, \
    const Matrix<T>& A, Matrix<T>& B, Int offset ); \
  template void AxpyTrapezoid \
  ( UpperOrLower uplo, S alpha, \
    const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B, Int offset ); \
  template void AxpyTrapezoid \
  ( UpperOrLower uplo, S alpha, \
    const SparseMatrix<T>& A, SparseMatrix<T>& B, Int offset ); \
  template void AxpyTrapezoid \
  ( UpperOrLower uplo, S alpha, \
    const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B, Int offset );

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T)

#define PROTO_COMPLEX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_TYPES(T,T)

#include "El/macros/Instantiate.h"

} // namespace El
