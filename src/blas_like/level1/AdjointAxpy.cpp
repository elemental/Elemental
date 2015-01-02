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
void AdjointAxpy( S alphaS, const Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("AdjointAxpy"))
    TransposeAxpy( alphaS, X, Y, true );
}

template<typename T,typename S>
void AdjointAxpy( S alphaS, const SparseMatrix<T>& X, SparseMatrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("AdjointAxpy"))
    TransposeAxpy( alphaS, X, Y, true );
}

template<typename T,typename S>
void AdjointAxpy
( S alphaS, const AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("AdjointAxpy"))
    TransposeAxpy( alphaS, X, Y, true );
}

template<typename T,typename S>
void AdjointAxpy
( S alphaS, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("AdjointAxpy"))
    TransposeAxpy( alphaS, X, Y, true );
}

#define PROTO_TYPES(T,S) \
  template void AdjointAxpy \
  ( S alpha, const Matrix<T>& A, Matrix<T>& B ); \
  template void AdjointAxpy \
  ( S alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void AdjointAxpy \
  ( S alpha, const SparseMatrix<T>& A, SparseMatrix<T>& B ); \
  template void AdjointAxpy \
  ( S alpha, const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B );

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
