/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_HPP
#define LAPACK_NORM_HPP

#include "elemental/lapack-like/Norm/EntrywiseOne.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"
#include "elemental/lapack-like/Norm/One.hpp"

#include "elemental/lapack-like/Norm/Nuclear.hpp"
#include "elemental/lapack-like/Norm/Two.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type
Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM )
{
#ifndef RELEASE
    PushCallStack("Norm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case ENTRYWISE_ONE_NORM:
        norm = EntrywiseOneNorm( A );
        break;
    case FROBENIUS_NORM: 
        norm = FrobeniusNorm( A );
        break;
    case INFINITY_NORM:
        norm = InfinityNorm( A );
        break;
    case MAX_NORM:
        norm = MaxNorm( A );
        break;
    case ONE_NORM:
        norm = OneNorm( A );
        break;
    // The following two norms make use of an SVD
    case NUCLEAR_NORM:
        norm = NuclearNorm( A );
        break;
    case TWO_NORM:
        norm = TwoNorm( A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
Norm( const DistMatrix<F,U,V>& A, NormType type=FROBENIUS_NORM )
{
#ifndef RELEASE
    PushCallStack("Norm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM: 
        norm = FrobeniusNorm( A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = EntrywiseOneNorm( A );
        break;
    case INFINITY_NORM:
        norm = InfinityNorm( A );
        break;
    case MAX_NORM:
        norm = MaxNorm( A );
        break;
    case ONE_NORM:
        norm = OneNorm( A );
        break;
    // The following norms make use of an SVD
    case NUCLEAR_NORM:
        norm = NuclearNorm( A );
        break;
    case TWO_NORM:
        norm = TwoNorm( A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_HPP
