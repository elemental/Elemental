/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_CONDITION_INFINITY_HPP
#define ELEM_LAPACK_CONDITION_INFINITY_HPP

#include "elemental/lapack-like/Inverse.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"

namespace elem {

template<typename F> 
inline BASE(F)
InfinityCondition( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("InfinityCondition");
#endif
    typedef BASE(F) Real;
    Matrix<F> B( A );
    const Real infNorm = InfinityNorm( B );
    try { Inverse( B ); }
    catch( SingularMatrixException& e ) 
    { return std::numeric_limits<Real>::infinity(); }
    const Real infNormInv = InfinityNorm( B );
    return infNorm*infNormInv;
}

template<typename F,Distribution U,Distribution V> 
inline BASE(F)
InfinityCondition( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("InfinityCondition");
#endif
    typedef BASE(F) Real;
    DistMatrix<F> B( A );
    const Real infNorm = InfinityNorm( B );
    try { Inverse( B ); }
    catch( SingularMatrixException& e ) 
    { return std::numeric_limits<Real>::infinity(); }
    const Real infNormInv = InfinityNorm( B );
    return infNorm*infNormInv;
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_CONDITION_INFINITY_HPP
