/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
Base<Field> Condition( const Matrix<Field>& A, NormType type )
{
    EL_DEBUG_CSE
    Base<Field> norm = 0;
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

template<typename Field>
Base<Field> Condition( const AbstractDistMatrix<Field>& A, NormType type )
{
    EL_DEBUG_CSE
    Base<Field> norm = 0;
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

#define PROTO(Field) \
  template Base<Field> Condition( const Matrix<Field>& A, NormType type ); \
  template Base<Field> \
  Condition( const AbstractDistMatrix<Field>& A, NormType type );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
