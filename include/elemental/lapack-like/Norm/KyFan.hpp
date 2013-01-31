/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_KYFAN_HPP
#define LAPACK_NORM_KYFAN_HPP

#include "elemental/lapack-like/SVD.hpp"

namespace elem {

template<typename F> 
inline typename Base<F>::type
KyFanNorm( const Matrix<F>& A, int k )
{
#ifndef RELEASE
    PushCallStack("KyFanNorm");
#endif
    if( k < 1 || k > std::min(A.Height(),A.Width()) )
        throw std::logic_error("Invalid index of KyFan norm");

    typedef typename Base<F>::type R;
    Matrix<F> B( A );
    Matrix<R> s;
    SingularValues( B, s );

    R norm = 0;
    for( int j=k-1; j>=0; --j )
        norm += s.Get(j,0);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
KyFanNorm( const DistMatrix<F,U,V>& A, int k )
{
#ifndef RELEASE
    PushCallStack("KyFanNorm");
#endif
    if( k < 1 || k > std::min(A.Height(),A.Width()) )
        throw std::logic_error("Invalid index of KyFan norm");

    typedef typename Base<F>::type R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SingularValues( B, s );

    R localNorm = 0;
    DistMatrix<R,VR,STAR> sTop( A.Grid() );
    LockedView( sTop, s, 0, 0, k, 1 );
    const int localHeight = sTop.LocalHeight();
    for( int j=localHeight-1; j>=0; --j )
        localNorm += sTop.GetLocal(j,0);
    R norm;
    mpi::AllReduce( &localNorm, &norm, 1, mpi::SUM, A.Grid().VRComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_KYFAN_HPP
