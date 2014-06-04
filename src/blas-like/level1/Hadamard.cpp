/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

// C(i,j) := A(i,j) B(i,j)

namespace El {

template<typename T> 
void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Hadamard"))
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    C.Resize( A.Height(), A.Width() );

    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            C.Set( i, j, A.Get(i,j)*B.Get(i,j) );
}

template<typename T,Dist U,Dist V> 
void Hadamard
( const DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B, DistMatrix<T,U,V>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Hadamard"))
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    if( A.Grid() != B.Grid() )
        LogicError("A and B must have the same grids");
    if( A.ColAlign() != B.ColAlign() || A.RowAlign() != B.RowAlign() )
        LogicError("A and B must be aligned");
    C.AlignWith( A );
    C.Resize( A.Height(), A.Width() );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const T alpha = A.GetLocal(iLoc,jLoc); 
            const T beta = B.GetLocal(iLoc,jLoc);
            C.SetLocal( iLoc, jLoc, alpha*beta );
        }
    }
}

template<typename T>
void Hadamard
( const AbstractDistMatrix<T>& A, 
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Dot"))
    if( A.DistData().colDist != B.DistData().colDist ||
        B.DistData().colDist != C.DistData().colDist ||
        A.DistData().rowDist != B.DistData().rowDist || 
        B.DistData().rowDist != C.DistData().rowDist )
        RuntimeError("{A,B,C} must have the same distribution");
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
        auto& BCast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(B); \
        auto& CCast = dynamic_cast<      DistMatrix<T,CDIST,RDIST>&>(C); \
        Hadamard( ACast, BCast, CCast );
    #include "El/core/GuardAndPayload.h"
}

#define DIST_PROTO(T,U,V) \
  template void Hadamard \
  ( const DistMatrix<T,U,V>& A, \
    const DistMatrix<T,U,V>& B, \
          DistMatrix<T,U,V>& C );

#define PROTO(T) \
  template void Hadamard \
  ( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C ); \
  template void Hadamard \
  ( const AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<T>& B, \
          AbstractDistMatrix<T>& C ); \
  DIST_PROTO(T,CIRC,CIRC) \
  DIST_PROTO(T,MC,  MR  ) \
  DIST_PROTO(T,MC,  STAR) \
  DIST_PROTO(T,MD,  STAR) \
  DIST_PROTO(T,MR,  MC  ) \
  DIST_PROTO(T,MR,  STAR) \
  DIST_PROTO(T,STAR,MC  ) \
  DIST_PROTO(T,STAR,MD  ) \
  DIST_PROTO(T,STAR,MR  ) \
  DIST_PROTO(T,STAR,STAR) \
  DIST_PROTO(T,STAR,VC  ) \
  DIST_PROTO(T,STAR,VR  ) \
  DIST_PROTO(T,VC,  STAR) \
  DIST_PROTO(T,VR,  STAR)

PROTO(Int)
PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
