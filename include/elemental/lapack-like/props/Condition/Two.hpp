/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONDITION_TWO_HPP
#define ELEM_CONDITION_TWO_HPP

#include ELEM_SVD_INC

namespace elem {

template<typename F> 
inline Base<F>
TwoCondition( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("TwoCondition"))
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

template<typename F,Dist U,Dist V> 
inline Base<F>
TwoCondition( const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("TwoCondition"))
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

#endif // ifndef ELEM_CONDITION_TWO_HPP
