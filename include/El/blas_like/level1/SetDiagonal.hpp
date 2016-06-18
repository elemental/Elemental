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
    DEBUG_CSE
    function<void(T&,T)> func
    ( []( T& beta, T gamma ) { beta = gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetRealPartOfDiagonal( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset )
{
    DEBUG_CSE
    function<void(T&,Base<T>)> func
    ( []( T& beta, Base<T> gamma ) { SetRealPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetImagPartOfDiagonal( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset )
{
    DEBUG_CSE
    function<void(T&,Base<T>)> func
    ( []( T& beta, Base<T> gamma ) { SetImagPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void SetDiagonal
( DistMatrix<T,U,V>& A, const ElementalMatrix<T>& d, Int offset )
{
    DEBUG_CSE
    function<void(T&,T)> func
    ( []( T& beta, T gamma ) { beta = gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void SetRealPartOfDiagonal
( DistMatrix<T,U,V>& A, const ElementalMatrix<Base<T>>& d, Int offset )
{
    DEBUG_CSE
    function<void(T&,Base<T>)> func
    ( []( T& beta, Base<T> gamma ) { SetRealPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void SetImagPartOfDiagonal
( DistMatrix<T,U,V>& A, const ElementalMatrix<Base<T>>& d, Int offset )
{
    DEBUG_CSE
    function<void(T&,Base<T>)> func
    ( []( T& beta, Base<T> gamma ) { SetImagPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetDiagonal
( ElementalMatrix<T>& A, const ElementalMatrix<T>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = static_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
      SetDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T>
void SetRealPartOfDiagonal
( ElementalMatrix<T>& A, const ElementalMatrix<Base<T>>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = static_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
      SetRealPartOfDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T>
void SetImagPartOfDiagonal
( ElementalMatrix<T>& A, const ElementalMatrix<Base<T>>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = static_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
      SetImagPartOfDiagonal( ACast, d, offset );
    #include <El/macros/GuardAndPayload.h>
}

} // namespace El

#endif // ifndef EL_BLAS_SETDIAGONAL_HPP
