/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_CONDITION_TWO_HPP
#define ELEM_LAPACK_CONDITION_TWO_HPP

#include "elemental/lapack-like/SVD.hpp"

namespace elem {

template<typename F> 
inline BASE(F)
TwoCondition( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("TwoCondition");
#endif
    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    SVD( B, s );

    R cond = 1;
    const Int numVals = s.Height();
    if( numVals > 0 )
        cond = s.Get(0,0) / s.Get(numVals-1,0);
    return cond;
}

template<typename F,Distribution U,Distribution V> 
inline BASE(F)
TwoCondition( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("TwoCondition");
#endif
    typedef Base<F> R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SVD( B, s );

    R cond = 1;
    const Int numVals = s.Height();
    if( numVals > 0 )
        cond = s.Get(0,0) / s.Get(numVals-1,0);
    return cond;
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_CONDITION_TWO_HPP
