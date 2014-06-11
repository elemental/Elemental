/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./Condition/Frobenius.hpp"
#include "./Condition/Infinity.hpp"
#include "./Condition/Max.hpp"
#include "./Condition/One.hpp"
#include "./Condition/Two.hpp"

namespace El {

template<typename F> 
Base<F> Condition( const Matrix<F>& A, NormType type )
{
    DEBUG_ONLY(CallStackEntry cse("Condition"))
    Base<F> norm = 0;
    switch( type )
    {
    case FROBENIUS_NORM:
        norm = FrobeniusCondition( A );
        break;
    case INFINITY_NORM:
        norm = InfinityCondition( A );
        break;
    case MAX_NORM:
        norm = MaxCondition( A );
        break;
    case ONE_NORM:
        norm = OneCondition( A );
        break;
    case TWO_NORM:
        norm = TwoCondition( A );
        break;
    default:
        LogicError("Invalid norm type for condition number");
    }
    return norm;
}

template<typename F,Dist U,Dist V> 
Base<F> Condition( const DistMatrix<F,U,V>& A, NormType type )
{
    DEBUG_ONLY(CallStackEntry cse("Condition"))
    Base<F> norm = 0;
    switch( type )
    {
    case FROBENIUS_NORM:
        norm = FrobeniusCondition( A );
        break;
    case INFINITY_NORM:
        norm = InfinityCondition( A );
        break;
    case MAX_NORM:
        norm = MaxCondition( A );
        break;
    case ONE_NORM:
        norm = OneCondition( A );
        break;
    case TWO_NORM:
        norm = TwoCondition( A );
        break;
    default:
        LogicError("Invalid norm type for condition number");
    }
    return norm;
}

#define PROTO_DIST(F,U,V) \
  template Base<F> Condition( const DistMatrix<F,U,V>& A, NormType type ); \
  template Base<F> FrobeniusCondition( const DistMatrix<F,U,V>& A ); \
  template Base<F> InfinityCondition( const DistMatrix<F,U,V>& A ); \
  template Base<F> MaxCondition( const DistMatrix<F,U,V>& A ); \
  template Base<F> OneCondition( const DistMatrix<F,U,V>& A ); \
  template Base<F> TwoCondition( const DistMatrix<F,U,V>& A );

#define PROTO(F) \
  template Base<F> Condition( const Matrix<F>& A, NormType type ); \
  template Base<F> FrobeniusCondition( const Matrix<F>& A ); \
  template Base<F> InfinityCondition( const Matrix<F>& A ); \
  template Base<F> MaxCondition( const Matrix<F>& A ); \
  template Base<F> OneCondition( const Matrix<F>& A ); \
  template Base<F> TwoCondition( const Matrix<F>& A ); \
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
