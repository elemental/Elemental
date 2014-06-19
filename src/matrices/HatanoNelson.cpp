/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Please see Section 36 of Trefethen and Embree's "Spectra and Pseudospectra"

template<typename F> 
void HatanoNelson
( Matrix<F>& A, Int n, F center, Base<F> radius, F g, bool periodic )
{
    DEBUG_ONLY(CallStackEntry cse("HatanoNelson"))
    if( n < 3 )
        LogicError("Hatano Nelson requires at least a 3x3 matrix");
    Zeros( A, n, n );
    auto d = A.GetDiagonal();
    MakeUniform( d, center, radius );
    A.SetDiagonal( d );
    SetDiagonal( A, Exp(g),   1 );
    SetDiagonal( A, Exp(-g), -1 );
    if( periodic )
    {
        A.Set( 0,   n-1, Exp(-g) );
        A.Set( n-1, 0,   Exp( g) );
    }
}

template<typename F,Dist U,Dist V>
void HatanoNelson
( DistMatrix<F,U,V>& A, Int n, F center, Base<F> radius, F g, bool periodic )
{
    DEBUG_ONLY(CallStackEntry cse("HatanoNelson"))
    if( n < 3 )
        LogicError("Hatano Nelson requires at least a 3x3 matrix");
    Zeros( A, n, n );
    auto d = A.GetDiagonal();
    MakeUniform( d, center, radius );
    A.SetDiagonal( d );
    SetDiagonal( A, Exp(g),   1 );
    SetDiagonal( A, Exp(-g), -1 );
    if( periodic )
    {
        A.Set( 0,   n-1, Exp(-g) );
        A.Set( n-1, 0,   Exp( g) );
    }
}

template<typename F,Dist U,Dist V>
void HatanoNelson
( BlockDistMatrix<F,U,V>& A, Int n, F center, Base<F> radius, F g, 
  bool periodic )
{
    DEBUG_ONLY(CallStackEntry cse("HatanoNelson"))
    if( n < 3 )
        LogicError("Hatano Nelson requires at least a 3x3 matrix");
    Zeros( A, n, n );
    auto d = A.GetDiagonal();
    MakeUniform( d, center, radius );
    A.SetDiagonal( d );
    SetDiagonal( A, Exp(g),   1 );
    SetDiagonal( A, Exp(-g), -1 );
    if( periodic )
    {
        A.Set( 0,   n-1, Exp(-g) );
        A.Set( n-1, 0,   Exp( g) );
    }
}

#define PROTO_DIST(F,U,V) \
  template void HatanoNelson \
  ( DistMatrix<F,U,V>& A, Int n, F center, Base<F> radius, F g, \
    bool periodic ); \
  template void HatanoNelson \
  ( BlockDistMatrix<F,U,V>& A, Int n, F center, Base<F> radius, F g, \
    bool periodic ); 

#define PROTO(F) \
  template void HatanoNelson \
  ( Matrix<F>& A, Int n, F center, Base<F> radius, F g, bool periodic ); \
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
