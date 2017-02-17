/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_GETDIAGONAL_HPP
#define EL_BLAS_GETDIAGONAL_HPP

namespace El {

template<typename T>
void GetDiagonal( const Matrix<T>& A, Matrix<T>& d, Int offset )
{
    EL_DEBUG_CSE
    function<T(const T&)> identity( []( const T& alpha ) { return alpha; } );
    GetMappedDiagonal( A, d, identity, offset );
}

template<typename T>
void GetRealPartOfDiagonal
( const Matrix<T>& A, Matrix<Base<T>>& d, Int offset )
{
    EL_DEBUG_CSE
    function<Base<T>(const T&)> realPart
    ( []( const T& alpha ) { return RealPart(alpha); } );
    GetMappedDiagonal( A, d, realPart, offset );
}

template<typename T>
void GetImagPartOfDiagonal
( const Matrix<T>& A, Matrix<Base<T>>& d, Int offset )
{
    EL_DEBUG_CSE
    function<Base<T>(const T&)> imagPart
    ( []( const T& alpha ) { return ImagPart(alpha); } );
    GetMappedDiagonal( A, d, imagPart, offset );
}

template<typename T>
Matrix<T> GetDiagonal( const Matrix<T>& A, Int offset )
{
    Matrix<T> d;
    GetDiagonal( A, d, offset );
    return d;
}

template<typename T>
Matrix<Base<T>> GetRealPartOfDiagonal( const Matrix<T>& A, Int offset )
{
    Matrix<Base<T>> d;
    GetRealPartOfDiagonal( A, d, offset );
    return d;
}

template<typename T>
Matrix<Base<T>> GetImagPartOfDiagonal( const Matrix<T>& A, Int offset )
{
    Matrix<Base<T>> d;
    GetImagPartOfDiagonal( A, d, offset );
    return d;
}

// TODO(poulson): SparseMatrix implementation

template<typename T,Dist U,Dist V,DistWrap wrap>
void GetDiagonal
( const DistMatrix<T,U,V,wrap>& A, AbstractDistMatrix<T>& d, Int offset )
{
    EL_DEBUG_CSE
    function<T(const T&)> identity( []( const T& alpha ) { return alpha; } );
    GetMappedDiagonal( A, d, identity, offset );
}

template<typename T,Dist U,Dist V,DistWrap wrap>
void GetRealPartOfDiagonal
( const DistMatrix<T,U,V,wrap>& A, AbstractDistMatrix<Base<T>>& d, Int offset )
{
    EL_DEBUG_CSE
    function<Base<T>(const T&)> realPart
    ( []( const T& alpha ) { return RealPart(alpha); } );
    GetMappedDiagonal( A, d, realPart, offset );
}

template<typename T,Dist U,Dist V,DistWrap wrap>
void GetImagPartOfDiagonal
( const DistMatrix<T,U,V,wrap>& A, AbstractDistMatrix<Base<T>>& d, Int offset )
{
    EL_DEBUG_CSE
    function<Base<T>(const T&)> imagPart
    ( []( const T& alpha ) { return ImagPart(alpha); } );
    GetMappedDiagonal( A, d, imagPart, offset );
}

template<typename T>
void GetDiagonal
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST,WRAP) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST && \
      A.Wrap() == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
      auto& ACast = static_cast<const DistMatrix<T,CDIST,RDIST,WRAP>&>(A); \
      GetDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T>
void GetRealPartOfDiagonal
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Base<T>>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST,WRAP) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST && \
      A.Wrap() == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
      auto& ACast = static_cast<const DistMatrix<T,CDIST,RDIST,WRAP>&>(A); \
      GetRealPartOfDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T>
void GetImagPartOfDiagonal
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Base<T>>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST,WRAP) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST && \
      A.Wrap() == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
      auto& ACast = static_cast<const DistMatrix<T,CDIST,RDIST,WRAP>&>(A); \
      GetImagPartOfDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T,Dist U,Dist V>
DistMatrix<T,DiagCol<U,V>(),DiagRow<U,V>()>
GetDiagonal( const DistMatrix<T,U,V>& A, Int offset )
{
    DistMatrix<T,DiagCol<U,V>(),DiagRow<U,V>()> d(A.Grid());
    GetDiagonal( A, d, offset );
    return d;
}

// TODO(poulson): Use a different output distribution after BlockMatrix has
// better operator= coverage
template<typename T,Dist U,Dist V>
DistMatrix<T,MC,STAR>
GetDiagonal( const DistMatrix<T,U,V,BLOCK>& A, Int offset )
{
    DistMatrix<T,MC,STAR> d(A.Grid());
    GetDiagonal( A, d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
DistMatrix<Base<T>,DiagCol<U,V>(),DiagRow<U,V>()>
GetRealPartOfDiagonal( const DistMatrix<T,U,V>& A, Int offset )
{
    DistMatrix<Base<T>,DiagCol<U,V>(),DiagRow<U,V>()> d(A.Grid());
    GetRealPartOfDiagonal( A, d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
DistMatrix<Base<T>,MC,STAR>
GetRealPartOfDiagonal( const DistMatrix<T,U,V,BLOCK>& A, Int offset )
{
    DistMatrix<Base<T>,MC,STAR> d(A.Grid());
    GetRealPartOfDiagonal( A, d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
DistMatrix<Base<T>,DiagCol<U,V>(),DiagRow<U,V>()>
GetImagPartOfDiagonal( const DistMatrix<T,U,V>& A, Int offset )
{
    DistMatrix<Base<T>,DiagCol<U,V>(),DiagRow<U,V>()> d(A.Grid());
    GetImagPartOfDiagonal( A, d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
DistMatrix<Base<T>,MC,STAR>
GetImagPartOfDiagonal( const DistMatrix<T,U,V,BLOCK>& A, Int offset )
{
    DistMatrix<Base<T>,MC,STAR> d(A.Grid());
    GetImagPartOfDiagonal( A, d, offset );
    return d;
}

// TODO: DistSparseMatrix implementation

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void GetDiagonal \
  ( const Matrix<T>& A, Matrix<T>& d, Int offset ); \
  EL_EXTERN template void GetRealPartOfDiagonal \
  ( const Matrix<T>& A, Matrix<Base<T>>& d, Int offset ); \
  EL_EXTERN template void GetImagPartOfDiagonal \
  ( const Matrix<T>& A, Matrix<Base<T>>& d, Int offset ); \
  EL_EXTERN template Matrix<T> GetDiagonal \
  ( const Matrix<T>& A, Int offset ); \
  EL_EXTERN template Matrix<Base<T>> GetRealPartOfDiagonal \
  ( const Matrix<T>& A, Int offset ); \
  EL_EXTERN template Matrix<Base<T>> GetImagPartOfDiagonal \
  ( const Matrix<T>& A, Int offset ); \
  EL_EXTERN template void GetDiagonal \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& d, \
    Int offset ); \
  EL_EXTERN template void GetRealPartOfDiagonal \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<Base<T>>& d, \
    Int offset ); \
  EL_EXTERN template void GetImagPartOfDiagonal \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<Base<T>>& d, \
    Int offset );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_GETDIAGONAL_HPP
