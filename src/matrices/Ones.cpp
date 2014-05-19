/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El/matrices/Ones.hpp"

namespace El {

template<typename T> 
void MakeOnes( Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, T(1) );
}

template<typename T,Dist U,Dist V>
void MakeOnes( DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    MakeOnes( A.Matrix() );
}

template<typename T,Dist U,Dist V>
void MakeOnes( BlockDistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    MakeOnes( A.Matrix() );
}

template<typename T>
void Ones( Matrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

template<typename T,Dist U,Dist V>
void Ones( DistMatrix<T,U,V>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

template<typename T,Dist U,Dist V>
void Ones( BlockDistMatrix<T,U,V>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

template<typename T>
void Ones( AbstractDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST 
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
      Ones( ACast, m, n ); 
    #include "El/core/GuardAndPayload.h"
}

template<typename T>
void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST 
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<BlockDistMatrix<T,CDIST,RDIST>&>(A); \
      Ones( ACast, m, n ); 
    #include "El/core/GuardAndPayload.h"
}

#define DISTPROTO(T,U,V) \
  template void Ones( DistMatrix<T,U,V>& A, Int m, Int n ); \
  template void Ones( BlockDistMatrix<T,U,V>& A, Int m, Int n );

#define PROTO(T) \
  template void Ones( Matrix<T>& A, Int m, Int n ); \
  DISTPROTO(T,CIRC,CIRC); \
  DISTPROTO(T,MC,  MR  ); \
  DISTPROTO(T,MC,  STAR); \
  DISTPROTO(T,MD,  STAR); \
  DISTPROTO(T,MR,  MC  ); \
  DISTPROTO(T,MR,  STAR); \
  DISTPROTO(T,STAR,MC  ); \
  DISTPROTO(T,STAR,MD  ); \
  DISTPROTO(T,STAR,MR  ); \
  DISTPROTO(T,STAR,STAR); \
  DISTPROTO(T,STAR,VC  ); \
  DISTPROTO(T,STAR,VR  ); \
  DISTPROTO(T,VC,  STAR); \
  DISTPROTO(T,VR,  STAR); \
  template void Ones( AbstractDistMatrix<T>& A, Int m, Int n ); \
  template void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n );

PROTO(Int);
#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
