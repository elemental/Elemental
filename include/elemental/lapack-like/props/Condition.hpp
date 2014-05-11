/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONDITION_HPP
#define ELEM_CONDITION_HPP

#include "./Condition/Frobenius.hpp"
#include "./Condition/Infinity.hpp"
#include "./Condition/Max.hpp"
#include "./Condition/One.hpp"
#include "./Condition/Two.hpp"

namespace elem {

template<typename F> 
inline Base<F>
Condition( const Matrix<F>& A, NormType type=TWO_NORM )
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
inline Base<F>
Condition( const DistMatrix<F,U,V>& A, NormType type=TWO_NORM )
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

} // namespace elem

#endif // ifndef ELEM_CONDITION_HPP
