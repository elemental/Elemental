/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_CONDITIONNUMBER_HPP
#define LAPACK_CONDITIONNUMBER_HPP

#include "elemental/lapack-like/SVD.hpp"

namespace elem {

template<typename F> 
inline typename Base<F>::type
ConditionNumber( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("ConditionNumber");
#endif
    typedef typename Base<F>::type R;

    Matrix<F> B( A );
    Matrix<R> s;
    SingularValues( B, s );

    R cond = 1;
    const int numVals = s.Height();
    if( numVals > 0 )
        cond = s.Get(0,0) / s.Get(numVals-1,0);
#ifndef RELEASE
    PopCallStack();
#endif
    return cond;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
ConditionNumber( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("ConditionNumber");
#endif
    typedef typename Base<F>::type R;

    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SingularValues( B, s );

    R cond = 1;
    const int numVals = s.Height();
    if( numVals > 0 )
        cond = s.Get(0,0) / s.Get(numVals-1,0);
#ifndef RELEASE
    PopCallStack();
#endif
    return cond;
}

} // namespace elem

#endif // ifndef LAPACK_CONDITIONNUMBER_HPP
