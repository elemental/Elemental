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
void SetDiagonal( Matrix<T>& A, const Matrix<T>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("SetDiagonal"))
    std::function<void(T&,T)> func
    ( []( T& beta, T gamma ) { beta = gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetRealPartOfDiagonal( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("SetRealPartOfDiagonal"))
    std::function<void(T&,Base<T>)> func
    ( []( T& beta, Base<T> gamma ) { SetRealPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetImagPartOfDiagonal( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("SetImagPartOfDiagonal"))
    std::function<void(T&,Base<T>)> func
    ( []( T& beta, Base<T> gamma ) { SetImagPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void SetDiagonal
( DistMatrix<T,U,V>& A, const AbstractDistMatrix<T>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("SetDiagonal"))
    std::function<void(T&,T)> func
    ( []( T& beta, T gamma ) { beta = gamma; } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void SetRealPartOfDiagonal
( DistMatrix<T,U,V>& A, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("SetRealPartOfDiagonal"))
    std::function<void(T&,Base<T>)> func
    ( []( T& beta, Base<T> gamma ) { SetRealPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T,Dist U,Dist V>
void SetImagPartOfDiagonal
( DistMatrix<T,U,V>& A, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("SetImagPartOfDiagonal"))
    std::function<void(T&,Base<T>)> func
    ( []( T& beta, Base<T> gamma ) { SetImagPart(beta,gamma); } );
    UpdateMappedDiagonal( A, d, func, offset );
}

template<typename T>
void SetDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
      SetDiagonal( ACast, d, offset );
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
void SetRealPartOfDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
      SetRealPartOfDiagonal( ACast, d, offset );
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
void SetImagPartOfDiagonal
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Base<T>>& d, Int offset )
{
    // Manual dynamic dispatch
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
      SetImagPartOfDiagonal( ACast, d, offset );
    #include "El/macros/GuardAndPayload.h"
}

#define PROTO_DIST(T,U,V) \
  template void SetDiagonal \
  ( DistMatrix<T,U,V>& A, const AbstractDistMatrix<T>& d, \
    Int offset ); \
  template void SetRealPartOfDiagonal \
  ( DistMatrix<T,U,V>& A, const AbstractDistMatrix<Base<T>>& d, \
    Int offset ); \
  template void SetImagPartOfDiagonal \
  ( DistMatrix<T,U,V>& A, const AbstractDistMatrix<Base<T>>& d, \
    Int offset );

#define PROTO(T) \
  template void SetDiagonal \
  ( Matrix<T>& A, const Matrix<T>& d, Int offset ); \
  template void SetRealPartOfDiagonal \
  ( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset ); \
  template void SetImagPartOfDiagonal \
  ( Matrix<T>& A, const Matrix<Base<T>>& d, Int offset ); \
  template void SetDiagonal \
  ( AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& d, \
    Int offset ); \
  template void SetRealPartOfDiagonal \
  ( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Base<T>>& d, \
    Int offset ); \
  template void SetImagPartOfDiagonal \
  ( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Base<T>>& d, \
    Int offset ); \
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

#include "El/macros/Instantiate.h"

} // namespace El
