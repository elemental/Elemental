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
void ApplyColPivots( Matrix<T>& A, const Matrix<Int>& pivots, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyColPivots");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
        if( pivots.Height() > A.Width() )
            LogicError("pivots cannot be longer than width of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = pivots.Height();
    for( Int j=0; j<numPivots; ++j )
    {
        const Int k = pivots.Get(j,0)-offset;
        T* Aj = A.Buffer(0,j);
        T* Ak = A.Buffer(0,k);
        for( Int i=0; i<height; ++i )
        {
            T temp = Aj[i];
            Aj[i] = Ak[i];
            Ak[i] = temp;
        }
    }
}

template<typename T>
void ApplyInverseColPivots
( Matrix<T>& A, const Matrix<Int>& pivots, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyInverseColPivots");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
        if( pivots.Height() > A.Width() )
            LogicError("pivots cannot be larger than width of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = pivots.Height();
    for( Int j=numPivots-1; j>=0; --j )
    {
        const Int k = pivots.Get(j,0)-offset;
        T* Aj = A.Buffer(0,j);
        T* Ak = A.Buffer(0,k);
        for( Int i=0; i<height; ++i )
        {
            T temp = Aj[i];
            Aj[i] = Ak[i];
            Ak[i] = temp;
        }
    }
}

template<typename T,Dist U,Dist V,Dist UPerm>
void ApplyColPivots
( DistMatrix<T,U,V>& A, const DistMatrix<Int,UPerm,STAR>& pivots, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyColPivots"))
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
    PermuteCols( A, perm, invPerm );
}

template<typename T,Dist U,Dist V,Dist UPerm>
void ApplyInverseColPivots
( DistMatrix<T,U,V>& A, const DistMatrix<Int,UPerm,STAR>& pivots, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyInverseColPivots"))
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
    PermuteCols( A, invPerm, perm );
}

#define PROTO_DIST(T,U,V) \
  template void ApplyColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MC,STAR>& pivots, Int offset ); \
  template void ApplyColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MD,STAR>& pivots, Int offset ); \
  template void ApplyColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MR,STAR>& pivots, Int offset ); \
  template void ApplyColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,VC,STAR>& pivots, Int offset ); \
  template void ApplyColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,VR,STAR>& pivots, Int offset ); \
  template void ApplyColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,STAR,STAR>& pivots, \
    Int offset ); \
  template void ApplyInverseColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MC,STAR>& pivots, Int offset ); \
  template void ApplyInverseColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MD,STAR>& pivots, Int offset ); \
  template void ApplyInverseColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,MR,STAR>& pivots, Int offset ); \
  template void ApplyInverseColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,VC,STAR>& pivots, Int offset ); \
  template void ApplyInverseColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,VR,STAR>& pivots, Int offset ); \
  template void ApplyInverseColPivots \
  ( DistMatrix<T,U,V>& A, const DistMatrix<Int,STAR,STAR>& pivots, \
    Int offset );

#define PROTO(T) \
  template void ApplyColPivots \
  ( Matrix<T>& A, const Matrix<Int>& pivots, Int offset ); \
  template void ApplyInverseColPivots \
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
