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
void ShiftDiagonal( Matrix<T>& A, S alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ShiftDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
    {
        const Int i = j-offset;
        if( i >= 0 && i < height )
            A.Update(i,j,alpha);
    }
}

template<typename T,typename S>
void ShiftDiagonal( AbstractDistMatrix<T>& A, S alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ShiftDiagonal"))
    const Int height = A.Height();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Int i = j-offset;
        if( i >= 0 && i < height && A.IsLocalRow(i) )
        {
            const Int iLoc = A.LocalRow(i);
            A.UpdateLocal( iLoc, jLoc, alpha );
        }
    }
}

template<typename T,typename S>
void ShiftDiagonal( AbstractBlockDistMatrix<T>& A, S alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ShiftDiagonal"))
    const Int height = A.Height();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Int i = j-offset;
        if( i >= 0 && i < height && A.IsLocalRow(i) )
        {
            const Int iLoc = A.LocalRow(i);
            A.UpdateLocal( iLoc, jLoc, alpha );
        }
    }
}

template<typename T,typename S>
void ShiftDiagonal( SparseMatrix<T>& A, S alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ShiftDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    A.Reserve( A.Capacity()+m );
    for( Int i=0; i<m; ++i )
    { 
        if( i+offset >= 0 && i+offset < n )
            A.QueueUpdate( i, i+offset, T(alpha) );
    }
    A.MakeConsistent();
}

template<typename T,typename S>
void ShiftDiagonal( DistSparseMatrix<T>& A, S alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ShiftDiagonal"))
    const Int mLocal = A.LocalHeight();
    const Int firstLocalRow = A.FirstLocalRow();
    const Int n = A.Width();
    A.Reserve( A.Capacity()+mLocal );
    for( Int iLocal=0; iLocal<mLocal; ++iLocal )
    {
        const Int i = iLocal+firstLocalRow;
        if( i+offset >= 0 && i+offset < n )
            A.QueueLocalUpdate( iLocal, i+offset, alpha );
    }
    A.MakeConsistent();
}

#define PROTO_TYPES(T,S) \
  template void ShiftDiagonal( Matrix<T>& A, S alpha, Int offset ); \
  template void ShiftDiagonal \
  ( AbstractDistMatrix<T>& A, S alpha, Int offset ); \
  template void ShiftDiagonal \
  ( AbstractBlockDistMatrix<T>& A, S alpha, Int offset ); \
  template void ShiftDiagonal( SparseMatrix<T>& A, S alpha, Int offset ); \
  template void ShiftDiagonal( DistSparseMatrix<T>& A, S alpha, Int offset ); 

#define PROTO_SAME(T) PROTO_TYPES(T,T) \

#define PROTO_INT(T) PROTO_SAME(T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_SAME(T)

#define PROTO_COMPLEX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_SAME(T)

#include "El/macros/Instantiate.h"

} // namespace El
