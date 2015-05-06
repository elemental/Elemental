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
void UpdateDiagonal( Matrix<T>& A, T alpha, const Matrix<T>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateDiagonal"))
    function<void(T&,T)> func
    ( [alpha]( T& beta, T gamma ) { beta += alpha*gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void UpdateRealPartOfDiagonal
( Matrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateRealPartOfDiagonal"))
    function<void(T&,Base<T>)> func
    ( [alpha]( T& beta, Base<T> gamma ) { UpdateRealPart(beta,alpha*gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void UpdateImagPartOfDiagonal
( Matrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateImagPartOfDiagonal"))
    function<void(T&,Base<T>)> func
    ( [alpha]( T& beta, Base<T> gamma ) { UpdateImagPart(beta,alpha*gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void UpdateDiagonal
( DistMatrix<T,U,V>& A, T alpha, const AbstractDistMatrix<T>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateDiagonal"))
    function<void(T&,T)> func
    ( [alpha]( T& beta, T gamma ) { beta += alpha*gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void UpdateRealPartOfDiagonal
( DistMatrix<T,U,V>& A, Base<T> alpha, const AbstractDistMatrix<Base<T>>& d, 
  Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateRealPartOfDiagonal"))
    function<void(T&,Base<T>)> func
    ( [alpha]( T& beta, Base<T> gamma ) { UpdateRealPart(beta,alpha*gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void UpdateImagPartOfDiagonal
( DistMatrix<T,U,V>& A, Base<T> alpha, const AbstractDistMatrix<Base<T>>& d, 
  Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateImagPartOfDiagonal"))
    function<void(T&,Base<T>)> func
    ( [alpha]( T& beta, Base<T> gamma ) { UpdateImagPart(beta,alpha*gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void UpdateDiagonal
( SparseMatrix<T>& A, T alpha, const Matrix<T>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateDiagonal"))
    function<void(T&,T)> func
    ( [alpha]( T& beta, T gamma ) { beta += alpha*gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void UpdateRealPartOfDiagonal
( SparseMatrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateRealPartOfDiagonal"))
    function<void(T&,Base<T>)> func
    ( [alpha]( T& beta, Base<T> gamma ) { UpdateRealPart(beta,alpha*gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void UpdateImagPartOfDiagonal
( SparseMatrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateImagPartOfDiagonal"))
    function<void(T&,Base<T>)> func
    ( [alpha]( T& beta, Base<T> gamma ) { UpdateImagPart(beta,alpha*gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void UpdateDiagonal
( DistSparseMatrix<T>& A, T alpha, const DistMultiVec<T>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateDiagonal"))
    function<void(T&,T)> func
    ( [alpha]( T& beta, T gamma ) { beta += alpha*gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void UpdateRealPartOfDiagonal
( DistSparseMatrix<T>& A, Base<T> alpha, 
  const DistMultiVec<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateRealPartOfDiagonal"))
    function<void(T&,Base<T>)> func
    ( [alpha]( T& beta, Base<T> gamma ) { UpdateRealPart(beta,alpha*gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void UpdateImagPartOfDiagonal
( DistSparseMatrix<T>& A, Base<T> alpha, 
  const DistMultiVec<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CSE cse("UpdateImagPartOfDiagonal"))
    function<void(T&,Base<T>)> func
    ( [alpha]( T& beta, Base<T> gamma ) { UpdateImagPart(beta,alpha*gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

#define PROTO_DIST(T,U,V) \
  template void UpdateDiagonal \
  ( DistMatrix<T,U,V>& A, T alpha, const AbstractDistMatrix<T>& d, \
    Int offset ); \
  template void UpdateRealPartOfDiagonal \
  ( DistMatrix<T,U,V>& A, Base<T> alpha, const AbstractDistMatrix<Base<T>>& d, \
    Int offset ); \
  template void UpdateImagPartOfDiagonal \
  ( DistMatrix<T,U,V>& A, Base<T> alpha, const AbstractDistMatrix<Base<T>>& d, \
    Int offset );

#define PROTO(T) \
  template void UpdateDiagonal \
  ( Matrix<T>& A, T alpha, const Matrix<T>& d, Int offset ); \
  template void UpdateRealPartOfDiagonal \
  ( Matrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset ); \
  template void UpdateImagPartOfDiagonal \
  ( Matrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset ); \
  template void UpdateDiagonal \
  ( SparseMatrix<T>& A, T alpha, const Matrix<T>& d, Int offset ); \
  template void UpdateRealPartOfDiagonal \
  ( SparseMatrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset ); \
  template void UpdateImagPartOfDiagonal \
  ( SparseMatrix<T>& A, Base<T> alpha, const Matrix<Base<T>>& d, Int offset ); \
  template void UpdateDiagonal \
  ( DistSparseMatrix<T>& A, T alpha, \
    const DistMultiVec<T>& d, Int offset ); \
  template void UpdateRealPartOfDiagonal \
  ( DistSparseMatrix<T>& A, Base<T> alpha, \
    const DistMultiVec<Base<T>>& d, Int offset ); \
  template void UpdateImagPartOfDiagonal \
  ( DistSparseMatrix<T>& A, Base<T> alpha, \
    const DistMultiVec<Base<T>>& d, Int offset ); \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,  MR  ) \
  PROTO_DIST(T,MC,  STAR) \
  PROTO_DIST(T,MD,  STAR) \
  PROTO_DIST(T,MR,  MC  ) \
  PROTO_DIST(T,MR,  STAR) \
  PROTO_DIST(T,STAR,MC  ) \
  PROTO_DIST(T,STAR,MD  ) \
  PROTO_DIST(T,STAR,MR  ) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC  ) \
  PROTO_DIST(T,STAR,VR  ) \
  PROTO_DIST(T,VC,  STAR) \
  PROTO_DIST(T,VR,  STAR)

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
