/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T>
void ApplySymmetricPivots
( UpperOrLower uplo, Matrix<T>& A, 
  const Matrix<Int>& p, bool conjugate, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplySymmetricPivots");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( p.Height() > A.Width() )
            LogicError("p cannot be longer than width of A");
        if( A.Height() != A.Width() )
            LogicError("A must be symmetric");
    )
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=0; k<n; ++k )
        SymmetricSwap( uplo, A, k, p.Get(k,0)-offset, conjugate );
}

template<typename T,Dist UPerm>
void ApplySymmetricPivots
( UpperOrLower uplo, DistMatrix<T>& A, 
  const DistMatrix<Int,UPerm,STAR>& pivots, bool conjugate, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplySymmetricPivots");
        if( pivots.Width() != 1 )
            LogicError("p must be a column vector");
        if( pivots.Height() > A.Width() )
            LogicError("p cannot be longer than width of A");
        if( A.Height() != A.Width() )
            LogicError("A must be symmetric");
    )
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=0; k<n; ++k )
        SymmetricSwap( uplo, A, k, pivots.Get(k,0)-offset, conjugate );
}

template<typename T>
void ApplyInverseSymmetricPivots
( UpperOrLower uplo, Matrix<T>& A, 
  const Matrix<Int>& p, bool conjugate, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyInverseSymmetricPivots");
        if( p.Width() != 1 )
            LogicError("pivots must be a column vector");
        if( p.Height() > A.Width() )
            LogicError("pivots cannot be longer than width of A");
        if( A.Height() != A.Width() )
            LogicError("A must be symmetric");
    )
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=n-1; k>=0; --k )
        SymmetricSwap( uplo, A, k, p.Get(k,0)-offset, conjugate );
}

template<typename T,Dist UPerm>
void ApplyInverseSymmetricPivots
( UpperOrLower uplo, DistMatrix<T>& A, 
  const DistMatrix<Int,UPerm,STAR>& pivots, bool conjugate, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyInverseSymmetricPivots");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
        if( pivots.Height() > A.Width() )
            LogicError("pivots cannot be longer than width of A");
        if( A.Height() != A.Width() )
            LogicError("A must be symmetric");
    )
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=n-1; k>=0; --k )
        SymmetricSwap( uplo, A, k, pivots.Get(k,0)-offset, conjugate );
}

#define PROTO_DIST(T,U) \
  template void ApplySymmetricPivots \
  ( UpperOrLower uplo, DistMatrix<T>& A, const DistMatrix<Int,U,STAR>& p, \
    bool conjugate, Int offset ); \
  template void ApplyInverseSymmetricPivots \
  ( UpperOrLower uplo, DistMatrix<T>& A, const DistMatrix<Int,U,STAR>& p, \
    bool conjugate, Int offset );

#define PROTO(T) \
  template void ApplySymmetricPivots \
  ( UpperOrLower uplo, Matrix<T>& A, const Matrix<Int>& p, \
    bool conjugate, Int offset ); \
  template void ApplyInverseSymmetricPivots \
  ( UpperOrLower uplo, Matrix<T>& A, const Matrix<Int>& p, \
    bool conjugate, Int offset ); \
  PROTO_DIST(T,MC  ) \
  PROTO_DIST(T,MD  ) \
  PROTO_DIST(T,MR  ) \
  PROTO_DIST(T,STAR) \
  PROTO_DIST(T,VC  ) \
  PROTO_DIST(T,VR  )

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
