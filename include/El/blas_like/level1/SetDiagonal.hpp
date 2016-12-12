/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_SETDIAGONAL_HPP
#define EL_BLAS_SETDIAGONAL_HPP

namespace El {

template<typename T>
void SetDiagonal( Matrix<T>& A, const Matrix<T>& d, Int offset )
{
    EL_DEBUG_CSE
    function<void(T&,const T&)> func
    ( []( T& beta, const T& gamma ) { beta = gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetRealPartOfDiagonal( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset )
{
    EL_DEBUG_CSE
    function<void(T&,const Base<T>&)> func
    ( []( T& beta, const Base<T>& gamma ) { SetRealPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetImagPartOfDiagonal( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset )
{
    EL_DEBUG_CSE
    function<void(T&,const Base<T>&)> func
    ( []( T& beta, const Base<T>& gamma ) { SetImagPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V,DistWrap wrap>
void SetDiagonal
( DistMatrix<T,U,V,wrap>& A, const AbstractDistMatrix<T>& d, Int offset )
{
    EL_DEBUG_CSE
    function<void(T&,const T&)> func
    ( []( T& beta, const T& gamma ) { beta = gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V,DistWrap wrap>
void SetRealPartOfDiagonal
( DistMatrix<T,U,V,wrap>& A, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    EL_DEBUG_CSE
    function<void(T&,const Base<T>&)> func
    ( []( T& beta, const Base<T>& gamma ) { SetRealPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V,DistWrap wrap>
void SetImagPartOfDiagonal
( DistMatrix<T,U,V,wrap>& A, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    EL_DEBUG_CSE
    function<void(T&,const Base<T>&)> func
    ( []( T& beta, const Base<T>& gamma ) { SetImagPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST,WRAP) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST && \
      A.Wrap() == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
      auto& ACast = static_cast<DistMatrix<T,CDIST,RDIST,WRAP>&>(A); \
      SetDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T>
void SetRealPartOfDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST,WRAP) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST && \
      A.Wrap() == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
      auto& ACast = static_cast<DistMatrix<T,CDIST,RDIST,WRAP>&>(A); \
      SetRealPartOfDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T>
void SetImagPartOfDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST,WRAP) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST && \
      A.Wrap() == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
      auto& ACast = static_cast<DistMatrix<T,CDIST,RDIST,WRAP>&>(A); \
      SetImagPartOfDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

} // namespace El

#endif // ifndef EL_BLAS_SETDIAGONAL_HPP
