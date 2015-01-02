/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void EntrywiseMap( Matrix<T>& A, std::function<T(T)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func(A.Get(i,j)) );
}

template<typename T>
void EntrywiseMap( SparseMatrix<T>& A, std::function<T(T)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    T* vBuf = A.ValueBuffer();
    const Int numEntries = A.NumEntries();
    for( Int k=0; k<numEntries; ++k )
        vBuf[k] = func(vBuf[k]);
}

template<typename T>
void EntrywiseMap( AbstractDistMatrix<T>& A, std::function<T(T)> func )
{ EntrywiseMap( A.Matrix(), func ); }

template<typename T>
void EntrywiseMap( AbstractBlockDistMatrix<T>& A, std::function<T(T)> func )
{ EntrywiseMap( A.Matrix(), func ); }

template<typename T>
void EntrywiseMap( DistSparseMatrix<T>& A, std::function<T(T)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    T* vBuf = A.ValueBuffer();
    const Int numLocalEntries = A.NumLocalEntries();
    for( Int k=0; k<numLocalEntries; ++k )
        vBuf[k] = func(vBuf[k]);
}

template<typename T>
void EntrywiseMap( DistMultiVec<T>& A, std::function<T(T)> func )
{ EntrywiseMap( A.Matrix(), func ); }

template<typename S,typename T>
void EntrywiseMap( const Matrix<S>& A, Matrix<T>& B, std::function<T(S)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            B.Set( i, j, func(A.Get(i,j)) );
}

template<typename S,typename T>
void EntrywiseMap
( const SparseMatrix<S>& A, SparseMatrix<T>& B, std::function<T(S)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    B.Empty();
    B.Resize( m, n );
    B.Reserve( numEntries );
    for( Int k=0; k<numEntries; ++k )
        B.QueueUpdate( A.Row(k), A.Col(k), func(A.Value(k)) );
    B.MakeConsistent();
}

template<typename S,typename T>
void EntrywiseMap
( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B, 
  std::function<T(S)> func )
{ 
    if( A.DistData().colDist == B.DistData().colDist &&
        A.DistData().rowDist == B.DistData().rowDist )
    {
        B.AlignWith( A.DistData() );
        B.Resize( A.Height(), A.Width() );
        EntrywiseMap( A.LockedMatrix(), B.Matrix(), func );
    }
    else
    {
        B.Resize( A.Height(), A.Width() );
        #define GUARD(CDIST,RDIST) \
          B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
        #define PAYLOAD(CDIST,RDIST) \
          DistMatrix<S,CDIST,RDIST> AProx(B.Grid()); \
          AProx.AlignWith( B.DistData() ); \
          Copy( A, AProx ); \
          EntrywiseMap( AProx.Matrix(), B.Matrix(), func );
        #include "El/macros/GuardAndPayload.h"
        #undef GUARD
        #undef PAYLOAD
    }
}

template<typename S,typename T>
void EntrywiseMap
( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B, 
  std::function<T(S)> func )
{ 
    if( A.DistData().colDist == B.DistData().colDist &&
        A.DistData().rowDist == B.DistData().rowDist )
    {
        B.AlignWith( A.DistData() );
        B.Resize( A.Height(), A.Width() );
        EntrywiseMap( A.LockedMatrix(), B.Matrix(), func );
    }
    else
    {
        B.Resize( A.Height(), A.Width() );
        #define GUARD(CDIST,RDIST) \
          B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
        #define PAYLOAD(CDIST,RDIST) \
          BlockDistMatrix<S,CDIST,RDIST> AProx(B.Grid()); \
          AProx.AlignWith( B.DistData() ); \
          Copy( A, AProx ); \
          EntrywiseMap( AProx.Matrix(), B.Matrix(), func );
        #include "El/macros/GuardAndPayload.h"
        #undef GUARD
        #undef PAYLOAD
    }
}

template<typename S,typename T>
void EntrywiseMap
( const DistSparseMatrix<S>& A, DistSparseMatrix<T>& B, 
  std::function<T(S)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numLocalEntries = A.NumLocalEntries();
    const Int firstLocalRow = A.FirstLocalRow();
    B.Empty();
    B.SetComm( A.Comm() );
    B.Resize( m, n );
    B.Reserve( numLocalEntries );
    for( Int k=0; k<numLocalEntries; ++k )
        B.QueueLocalUpdate
        ( A.Row(k)-firstLocalRow, A.Col(k), func(A.Value(k)) );
    B.MakeConsistent();
}

template<typename S,typename T>
void EntrywiseMap
( const DistMultiVec<S>& A, DistMultiVec<T>& B, std::function<T(S)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.SetComm( A.Comm() );
    B.Resize( m, n );
    EntrywiseMap( A.LockedMatrix(), B.Matrix(), func );
}

#define PROTO_TYPES(S,T) \
  template void EntrywiseMap \
  ( const Matrix<S>& A, Matrix<T>& B, std::function<T(S)> func ); \
  template void EntrywiseMap \
  ( const SparseMatrix<S>& A, SparseMatrix<T>& B, std::function<T(S)> func ); \
  template void EntrywiseMap \
  ( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B, \
    std::function<T(S)> func ); \
  template void EntrywiseMap \
  ( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B, \
    std::function<T(S)> func ); \
  template void EntrywiseMap \
  ( const DistSparseMatrix<S>& A, DistSparseMatrix<T>& B, \
    std::function<T(S)> func ); \
  template void EntrywiseMap \
  ( const DistMultiVec<S>& A, DistMultiVec<T>& B, std::function<T(S)> func );

#define PROTO(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,float) \
  PROTO_TYPES(T,double) \
  PROTO_TYPES(T,Complex<float>) \
  PROTO_TYPES(T,Complex<double>) \
  template void EntrywiseMap( Matrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap( SparseMatrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap \
  ( AbstractDistMatrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap \
  ( AbstractBlockDistMatrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap \
  ( DistSparseMatrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap \
  ( DistMultiVec<T>& A, std::function<T(T)> func );

#include "El/macros/Instantiate.h"

} // namespace El
