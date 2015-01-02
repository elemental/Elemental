/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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

template<typename F> 
Base<F> Condition( const AbstractDistMatrix<F>& A, NormType type )
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

#define PROTO(F) \
  template Base<F> Condition( const Matrix<F>& A, NormType type ); \
  template Base<F> Condition( const AbstractDistMatrix<F>& A, NormType type ); \
  template Base<F> FrobeniusCondition( const Matrix<F>& A ); \
  template Base<F> FrobeniusCondition( const AbstractDistMatrix<F>& A ); \
  template Base<F> InfinityCondition( const Matrix<F>& A ); \
  template Base<F> InfinityCondition( const AbstractDistMatrix<F>& A ); \
  template Base<F> MaxCondition( const Matrix<F>& A ); \
  template Base<F> MaxCondition( const AbstractDistMatrix<F>& A ); \
  template Base<F> OneCondition( const Matrix<F>& A ); \
  template Base<F> OneCondition( const AbstractDistMatrix<F>& A ); \
  template Base<F> TwoCondition( const Matrix<F>& A ); \
  template Base<F> TwoCondition( const AbstractDistMatrix<F>& A );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
