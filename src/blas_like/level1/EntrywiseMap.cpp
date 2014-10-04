/*
   Copyright (c) 2009-2014, Jack Poulson
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
void EntrywiseMap( AbstractDistMatrix<T>& A, std::function<T(T)> func )
{ EntrywiseMap( A.Matrix(), func ); }

template<typename T>
void EntrywiseMap( AbstractBlockDistMatrix<T>& A, std::function<T(T)> func )
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
( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B, 
  std::function<T(S)> func )
{ 
    B.AlignWith( A.DistData() );
    B.Resize( A.Height(), A.Width() );
    EntrywiseMap( A.LockedMatrix(), B.Matrix(), func ); 
}

template<typename S,typename T>
void EntrywiseMap
( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B, 
  std::function<T(S)> func )
{ 
    B.AlignWith( A.DistData() );
    B.Resize( A.Height(), A.Width() );
    EntrywiseMap( A.LockedMatrix(), B.Matrix(), func ); 
}

#define PROTO(T) \
  template void EntrywiseMap( Matrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap \
  ( AbstractDistMatrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap \
  ( AbstractBlockDistMatrix<T>& A, std::function<T(T)> func );

#define PROTO_TYPES(S,T) \
  template void EntrywiseMap \
  ( const Matrix<S>& A, Matrix<T>& B, std::function<T(S)> func ); \
  template void EntrywiseMap \
  ( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B, \
    std::function<T(S)> func ); \
  template void EntrywiseMap \
  ( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B, \
    std::function<T(S)> func );

#define PROTO_INT(T) \
  PROTO(T) \
  PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO(T) \
  PROTO_TYPES(Int,T) \
  PROTO_TYPES(T,T)

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  PROTO_TYPES(float,double) \
  PROTO_TYPES(float,Complex<double>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  PROTO_TYPES(double,float) \
  PROTO_TYPES(double,Complex<float>)

#define PROTO_COMPLEX(T) \
  PROTO(T) \
  PROTO_TYPES(Int,T) \
  PROTO_TYPES(Base<T>,T) 

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  PROTO_TYPES(Complex<float>,Complex<double>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  PROTO_TYPES(Complex<double>,Complex<float>)

#include "El/macros/Instantiate.h"

} // namespace El
