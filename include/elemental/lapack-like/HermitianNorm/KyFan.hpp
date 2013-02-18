/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANNORM_KYFAN_HPP
#define LAPACK_HERMITIANNORM_KYFAN_HPP

#include "elemental/blas-like/level1/MakeHermitian.hpp"

#ifndef WITHOUT_PMRRR
#include "elemental/lapack-like/HermitianSVD.hpp"
#endif // ifndef WITHOUT_PMRRR
#include "elemental/lapack-like/SVD.hpp"

namespace elem {

template<typename F> 
inline typename Base<F>::type
HermitianKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, int k )
{
#ifndef RELEASE
    PushCallStack("HermitianKyFanNorm");
#endif
    if( k < 1 || k > std::min(A.Height(),A.Width()) )
        throw std::logic_error("Invalid index of KyFan norm");

    typedef typename Base<F>::type R;
    Matrix<F> B( A );
    Matrix<R> s;
// TODO: Enable sequential MRRR
/*
#ifndef WITHOUT_PMRRR
    HermitianSingularValues( uplo, B, s );
#else
    MakeHermitian( uplo, B );
    SingularValues( B, s );
#endif // ifndef WITHOUT_PMRRR
*/
    MakeHermitian( uplo, B );
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
HermitianKyFanNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A, int k )
{
#ifndef RELEASE
    PushCallStack("HermitianKyFanNorm");
#endif
    if( k < 1 || k > std::min(A.Height(),A.Width()) )
        throw std::logic_error("Invalid index of KyFan norm");

    typedef typename Base<F>::type R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
#ifndef WITHOUT_PMRRR
    HermitianSingularValues( uplo, B, s );
#else
    MakeHermitian( uplo, B );
    SingularValues( B, s );
#endif // ifndef WITHOUT_PMRRR

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

#endif // ifndef LAPACK_HERMITIANNORM_KYFAN_HPP
