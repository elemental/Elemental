/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T,typename S>
void Axpy( S alphaS, const Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("Axpy"))
    const T alpha = T(alphaS);
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( X.Height()==1 || X.Width()==1 )
    {
        const Int XLength = ( X.Width()==1 ? X.Height() : X.Width() );
        const Int XStride = ( X.Width()==1 ? 1          : X.LDim() );
        const Int YStride = ( Y.Width()==1 ? 1          : Y.LDim() );
        DEBUG_ONLY(
            const Int YLength = ( Y.Width()==1 ? Y.Height() : Y.Width() );
            if( XLength != YLength )
                LogicError("Nonconformal Axpy");
        )
        blas::Axpy
        ( XLength, alpha, X.LockedBuffer(), XStride, Y.Buffer(), YStride );
    }
    else
    {
        DEBUG_ONLY(
            if( X.Height() != Y.Height() || X.Width() != Y.Width() )
                LogicError("Nonconformal Axpy");
        )
        if( X.Width() <= X.Height() )
        {
            for( Int j=0; j<X.Width(); ++j )
            {
                blas::Axpy
                ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( Int i=0; i<X.Height(); ++i )
            {
                blas::Axpy
                ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                    Y.Buffer(i,0),       Y.LDim() );
            }
        }
    }
}

template<typename T,typename S>
void Axpy( S alphaS, const SparseMatrix<T>& X, SparseMatrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("Axpy"))
    if( X.Height() != Y.Height() || X.Width() != Y.Width() )
        LogicError("X and Y must have the same dimensions");
    const T alpha = T(alphaS);
    const Int numEntries = X.NumEntries();
    Y.Reserve( Y.NumEntries()+numEntries );
    for( Int k=0; k<numEntries; ++k ) 
        Y.QueueUpdate( X.Row(k), X.Col(k), alpha*X.Value(k) );
    Y.MakeConsistent();
}

template<typename T,typename S>
void Axpy( S alphaS, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("Axpy");
        AssertSameGrids( X, Y );
    )
    const T alpha = T(alphaS);

    const DistData XDistData = X.DistData();
    const DistData YDistData = Y.DistData();

    if( XDistData == YDistData )
    {
        Axpy( alpha, X.LockedMatrix(), Y.Matrix() );
    }
    else
    {
        std::unique_ptr<AbstractDistMatrix<T>> 
          XCopy( Y.Construct(Y.Grid(),Y.Root()) );
        XCopy->AlignWith( YDistData );
        Copy( X, *XCopy );
        Axpy( alpha, XCopy->LockedMatrix(), Y.Matrix() );
    }
}

template<typename T,typename S>
void Axpy( S alphaS, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("Axpy"))
    if( X.Height() != Y.Height() || X.Width() != Y.Width() )
        LogicError("X and Y must have the same dimensions");
    if( X.Comm() != Y.Comm() )
        LogicError("X and Y must have the same communicator");
    const T alpha = T(alphaS);
    const Int numLocalEntries = X.NumLocalEntries();
    const Int firstLocalRow = X.FirstLocalRow();
    Y.Reserve( Y.NumLocalEntries()+numLocalEntries );
    for( Int k=0; k<numLocalEntries; ++k ) 
        Y.QueueLocalUpdate
        ( X.Row(k)-firstLocalRow, X.Col(k), alpha*X.Value(k) );
    Y.MakeConsistent();
}

template<typename T,typename S>
void Axpy( S alpha, const DistMultiVec<T>& X, DistMultiVec<T>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("Axpy");
        if( !mpi::Congruent( X.Comm(), Y.Comm() ) )
            LogicError("X and Y must have congruent communicators");
        if( X.Height() != Y.Height() )
            LogicError("X and Y must be the same height");
        if( X.Width() != Y.Width() )
            LogicError("X and Y must be the same width");
    )
    const int localHeight = X.LocalHeight();
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            Y.UpdateLocal( iLocal, j, T(alpha)*X.GetLocal(iLocal,j) );
}

#define PROTO_TYPES(T,S) \
  template void Axpy( S alpha, const Matrix<T>& A, Matrix<T>& B ); \
  template void Axpy \
  ( S alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void Axpy \
  ( S alpha, const SparseMatrix<T>& A, SparseMatrix<T>& B ); \
  template void Axpy \
  ( S alpha, const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B ); \
  template void Axpy( S alpha, const DistMultiVec<T>& X, DistMultiVec<T>& Y );

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
