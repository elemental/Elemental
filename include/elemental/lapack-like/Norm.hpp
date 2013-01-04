/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Norm/Util.hpp"

#include "./Norm/One.hpp"
#include "./Norm/Infinity.hpp"
#include "./Norm/Max.hpp"

#include "./Norm/Nuclear.hpp"
#include "./Norm/Frobenius.hpp"
#include "./Norm/Two.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type
Norm( const Matrix<F>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("Norm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    case ONE_NORM:
        norm = internal::OneNorm( A );
        break;
    case INFINITY_NORM:
        norm = internal::InfinityNorm( A );
        break;
    case MAX_NORM:
        norm = internal::MaxNorm( A );
        break;
    case NUCLEAR_NORM:
        norm = internal::NuclearNorm( A );
        break;
    case FROBENIUS_NORM: 
        norm = internal::FrobeniusNorm( A );
        break;
    case TWO_NORM:
        norm = internal::TwoNorm( A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
Norm( const DistMatrix<F,U,V>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("Norm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    case ONE_NORM:
        norm = internal::OneNorm( A );
        break;
    case INFINITY_NORM:
        norm = internal::InfinityNorm( A );
        break;
    case MAX_NORM:
        norm = internal::MaxNorm( A );
        break;
    case NUCLEAR_NORM:
        norm = internal::NuclearNorm( A );
        break;
    case FROBENIUS_NORM: 
        norm = internal::FrobeniusNorm( A );
        break;
    case TWO_NORM:
        norm = internal::TwoNorm( A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem
