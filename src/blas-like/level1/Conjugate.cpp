/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename Real>
void Conjugate( Matrix<Real>& A )
{ }

template<typename Real>
void Conjugate( Matrix<Complex<Real>>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Conjugate (in-place)"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set(i,j,Conj(A.Get(i,j)));
}

template<typename T>
void Conjugate( const Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Conjugate"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            B.Set(i,j,Conj(A.Get(i,j)));
}

template<typename T>
void Conjugate( AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Conjugate (in-place)"))
    Conjugate( A.Matrix() );
}

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Conjugate"))
    B = A;
    Conjugate( B );
}

template<typename T>
void Conjugate( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Axpy"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(B); \
        Conjugate( ACast, BCast );
    #include "El/core/NestedGuardAndPayload.h"
}

#define DIST_PROTO_INNER(T,U,V,W,Z) \
  template void Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B ); \

#define DIST_PROTO(T,U,V) \
  DIST_PROTO_INNER(T,U,V,CIRC,CIRC); \
  DIST_PROTO_INNER(T,U,V,MC,  MR  ); \
  DIST_PROTO_INNER(T,U,V,MC,  STAR); \
  DIST_PROTO_INNER(T,U,V,MD,  STAR); \
  DIST_PROTO_INNER(T,U,V,MR,  MC  ); \
  DIST_PROTO_INNER(T,U,V,MR,  STAR); \
  DIST_PROTO_INNER(T,U,V,STAR,MC  ); \
  DIST_PROTO_INNER(T,U,V,STAR,MD  ); \
  DIST_PROTO_INNER(T,U,V,STAR,MR  ); \
  DIST_PROTO_INNER(T,U,V,STAR,STAR); \
  DIST_PROTO_INNER(T,U,V,STAR,VC  ); \
  DIST_PROTO_INNER(T,U,V,STAR,VR  ); \
  DIST_PROTO_INNER(T,U,V,VC,  STAR); \
  DIST_PROTO_INNER(T,U,V,VR,  STAR); \

#define PROTO(T) \
  template void Conjugate( Matrix<T>& A ); \
  template void Conjugate( const Matrix<T>& A, Matrix<T>& B ); \
  template void Conjugate( AbstractDistMatrix<T>& A ); \
  template void Conjugate \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  DIST_PROTO(T,CIRC,CIRC); \
  DIST_PROTO(T,MC,  MR  ); \
  DIST_PROTO(T,MC,  STAR); \
  DIST_PROTO(T,MD,  STAR); \
  DIST_PROTO(T,MR,  MC  ); \
  DIST_PROTO(T,MR,  STAR); \
  DIST_PROTO(T,STAR,MC  ); \
  DIST_PROTO(T,STAR,MD  ); \
  DIST_PROTO(T,STAR,MR  ); \
  DIST_PROTO(T,STAR,STAR); \
  DIST_PROTO(T,STAR,VC  ); \
  DIST_PROTO(T,STAR,VR  ); \
  DIST_PROTO(T,VC  ,STAR); \
  DIST_PROTO(T,VR  ,STAR);

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
