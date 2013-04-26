/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_TWOUPPERBOUND_HPP
#define LAPACK_NORM_TWOUPPERBOUND_HPP

#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"
#include "elemental/lapack-like/Norm/One.hpp"

namespace elem {

template<typename F>
inline BASE(F)
TwoNormUpperBound( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("TwoNormUpperBound");
#endif
    typedef BASE(F) R;
    const R m = A.Height();
    const R n = A.Width();

    const R maxNorm = MaxNorm( A );
    const R oneNorm = OneNorm( A );
    const R infNorm = InfinityNorm( A );

    R upperBound = std::min( Sqrt(m*n)*maxNorm, Sqrt(m)*infNorm );
    upperBound = std::min( upperBound, Sqrt(n)*oneNorm );
    upperBound = std::min( upperBound, Sqrt( oneNorm*infNorm ) );
    return upperBound;
}

template<typename F> 
inline BASE(F)
TwoNormUpperBound( const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("TwoNormUpperBound");
#endif
    typedef BASE(F) R;
    const R m = A.Height();
    const R n = A.Width();

    const R maxNorm = MaxNorm( A );
    const R oneNorm = OneNorm( A );
    const R infNorm = InfinityNorm( A );

    R upperBound = std::min( Sqrt(m*n)*maxNorm, Sqrt(m)*infNorm );
    upperBound = std::min( upperBound, Sqrt(n)*oneNorm );
    upperBound = std::min( upperBound, Sqrt( oneNorm*infNorm ) );
    return upperBound;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_TWOUPPERBOUND_HPP
