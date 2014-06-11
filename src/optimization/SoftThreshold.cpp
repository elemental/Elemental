/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F>
F SoftThreshold( F alpha, Base<F> tau )
{
    DEBUG_ONLY(
        CallStackEntry cse("SoftThreshold");
        if( tau < 0 )
            LogicError("Negative threshold does not make sense");
    )
    const Base<F> scale = Abs(alpha);
    return ( scale <= tau ? F(0) : alpha-(alpha/scale)*tau );
}

template<typename F>
void SoftThreshold( Matrix<F>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("SoftThreshold"))
    if( relative )
        tau *= MaxNorm(A);
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            A.Set( i, j, SoftThreshold(A.Get(i,j),tau) );
}

template<typename F,Dist U,Dist V>
void SoftThreshold( DistMatrix<F,U,V>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CallStackEntry cse("SoftThreshold"))
    if( relative )
        tau *= MaxNorm(A);
    SoftThreshold( A.Matrix(), tau, false );
}

#define PROTO_DIST(F,U,V) \
  template void SoftThreshold \
  ( DistMatrix<F,U,V>& A, Base<F> tau, bool relative );

#define PROTO(F) \
  template F SoftThreshold( F alpha, Base<F> tau ); \
  template void SoftThreshold( Matrix<F>& A, Base<F> tau, bool relative ); \
  PROTO_DIST(F,CIRC,CIRC) \
  PROTO_DIST(F,MC,  MR  ) \
  PROTO_DIST(F,MC,  STAR) \
  PROTO_DIST(F,MD,  STAR) \
  PROTO_DIST(F,MR,  MC  ) \
  PROTO_DIST(F,MR,  STAR) \
  PROTO_DIST(F,STAR,MC  ) \
  PROTO_DIST(F,STAR,MD  ) \
  PROTO_DIST(F,STAR,MR  ) \
  PROTO_DIST(F,STAR,STAR) \
  PROTO_DIST(F,STAR,VC  ) \
  PROTO_DIST(F,STAR,VR  ) \
  PROTO_DIST(F,VC,  STAR) \
  PROTO_DIST(F,VR,  STAR)

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
