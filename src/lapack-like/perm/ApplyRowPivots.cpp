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
void ApplyRowPivots( Matrix<T>& A, const Matrix<Int>& pivots, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyRowPivots");
        if( pivots.Width() != 1 )
            LogicError("p must be a column vector");
        if( pivots.Height() > A.Height() )
            LogicError("p cannot be larger than height of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = pivots.Height();
    const Int ldim = A.LDim();
    for( Int i=0; i<numPivots; ++i )
    {
        const Int k = pivots.Get(i,0)-offset;
        T* Ai = A.Buffer(i,0);
        T* Ak = A.Buffer(k,0);
        for( Int j=0; j<width; ++j )
        {
            T temp = Ai[j*ldim];
            Ai[j*ldim] = Ak[j*ldim];
            Ak[j*ldim] = temp;
        }
    }
}

template<typename T>
void ApplyInverseRowPivots
( Matrix<T>& A, const Matrix<Int>& pivots, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyInverseRowPivots");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
        if( pivots.Height() > A.Height() )
            LogicError("pivots cannot be larger than height of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = pivots.Height();
    const Int ldim = A.LDim();
    for( Int i=numPivots-1; i>=0; --i )
    {
        const Int k = pivots.Get(i,0)-offset;
        T* Ai = A.Buffer(i,0);
        T* Ak = A.Buffer(k,0);
        for( Int j=0; j<width; ++j )
        {
            T temp = Ai[j*ldim];
            Ai[j*ldim] = Ak[j*ldim];
            Ak[j*ldim] = temp;
        }
    }
}

template<typename T,Dist U,Dist V,Dist UPerm>
void ApplyRowPivots
( DistMatrix<T,U,V>& A, const DistMatrix<Int,UPerm,STAR>& pivots, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyRowPivots"))
    DistMatrix<Int,UPerm,STAR> perm(pivots.Grid()),
                               invPerm(pivots.Grid());
    if( pivots.Height() == A.Width() )
    {
        PivotsToInversePermutation( pivots, invPerm, offset );
        InvertPermutation( invPerm, perm );
    }
    else
    {
        PivotsToPartialPermutation( pivots, perm, invPerm, offset );
    }
    PermuteRows( A, perm, invPerm );
}

template<typename T,Dist U,Dist V,Dist UPerm>
void ApplyInverseRowPivots
( DistMatrix<T,U,V>& A, const DistMatrix<Int,UPerm,STAR>& pivots, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyInverseRowPivots"))
    DistMatrix<Int,UPerm,STAR> perm(pivots.Grid()),
                               invPerm(pivots.Grid());
    if( pivots.Height() == A.Width() )
    {
        PivotsToInversePermutation( pivots, invPerm, offset );
        InvertPermutation( invPerm, perm );
    }
    else
    {
        PivotsToPartialPermutation( pivots, perm, invPerm, offset );
    }
    PermuteRows( A, invPerm, perm );
}

#define PROTO_DIST(T,U,V) \
  template void ApplyRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MC,STAR>& pivots, Int offset ); \
  template void ApplyRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MD,STAR>& pivots, Int offset ); \
  template void ApplyRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MR,STAR>& pivots, Int offset ); \
  template void ApplyRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,VC,STAR>& pivots, Int offset ); \
  template void ApplyRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,VR,STAR>& pivots, Int offset ); \
  template void ApplyRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,STAR,STAR>& pivots, \
    Int offset ); \
  template void ApplyInverseRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MC,STAR>& pivots, Int offset ); \
  template void ApplyInverseRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MD,STAR>& pivots, Int offset ); \
  template void ApplyInverseRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MR,STAR>& pivots, Int offset ); \
  template void ApplyInverseRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,VC,STAR>& pivots, Int offset ); \
  template void ApplyInverseRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,VR,STAR>& pivots, Int offset ); \
  template void ApplyInverseRowPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,STAR,STAR>& pivots, \
    Int offset );

#define PROTO(T) \
  template void ApplyRowPivots \
  ( Matrix<T>& A, const Matrix<Int>& pivots, Int offset ); \
  template void ApplyInverseRowPivots \
  ( Matrix<T>& A, const Matrix<Int>& pivots, Int offset ); \
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

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
