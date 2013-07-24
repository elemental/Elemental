/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_NORM_TWOLOWERBOUND_HPP
#define ELEM_LAPACK_NORM_TWOLOWERBOUND_HPP

#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"
#include "elemental/lapack-like/Norm/One.hpp"

namespace elem {

template<typename F>
inline BASE(F)
TwoNormLowerBound( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("TwoNormLowerBound");
#endif
    typedef BASE(F) R;
    const R m = A.Height();
    const R n = A.Width();

    const R maxNorm = MaxNorm( A );
    const R oneNorm = OneNorm( A );
    const R infNorm = InfinityNorm( A );
    const R frobNorm = FrobeniusNorm( A );
    R lowerBound = std::max( maxNorm, infNorm/Sqrt(n) );
    lowerBound = std::max( lowerBound, oneNorm/Sqrt(m) );
    lowerBound = std::max( lowerBound, frobNorm/Sqrt(std::min(m,n)) );
    return lowerBound;
}

template<typename F> 
inline BASE(F)
TwoNormLowerBound( const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("TwoNormLowerBound");
#endif
    typedef BASE(F) R;
    const R m = A.Height();
    const R n = A.Width();

    const R maxNorm = MaxNorm( A );
    const R oneNorm = OneNorm( A );
    const R infNorm = InfinityNorm( A );
    const R frobNorm = FrobeniusNorm( A );
    R lowerBound = std::max( maxNorm, infNorm/Sqrt(n) );
    lowerBound = std::max( lowerBound, oneNorm/Sqrt(m) );
    lowerBound = std::max( lowerBound, frobNorm/Sqrt(std::min(m,n)) );
    return lowerBound;
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_NORM_TWOLOWERBOUND_HPP
