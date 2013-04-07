/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_SCHATTEN_HPP
#define LAPACK_NORM_SCHATTEN_HPP

#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/blas-like/level1/MakeSymmetric.hpp"
#include "elemental/lapack-like/SVD.hpp"

namespace elem {

template<typename F> 
inline typename Base<F>::type
SchattenNorm( const Matrix<F>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("SchattenNorm");
#endif
    typedef typename Base<F>::type R;
    Matrix<F> B( A );
    Matrix<R> s;
    SVD( B, s );

    // TODO: Think of how to make this more stable
    const int k = s.Height();
    R sum = 0;
    for( int j=k-1; j>=0; --j )
        sum += Pow( s.Get(j,0), p ); 
    const R norm = Pow( sum, 1/p ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
HermitianSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("HermitianSchattenNorm");
#endif
    typedef typename Base<F>::type R;

    Matrix<F> B( A );
    Matrix<R> s;
    HermitianSVD( uplo, B, s );

    // TODO: Think of how to make this more stable
    const int k = s.Height();
    R sum = 0;
    for( int j=k-1; j>=0; --j )
        sum += Pow( s.Get(j,0), p );
    const R norm = Pow( sum, 1/p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
SymmetricSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("SymmetricSchattenNorm");
#endif
    typedef typename Base<F>::type R;
    Matrix<F> B( A );
    Matrix<R> s;
    MakeSymmetric( uplo, B );
    SVD( B, s );

    // TODO: Think of how to make this more stable
    const int k = s.Height();
    R sum = 0;
    for( int j=0; j<k; ++j )
        sum += Pow( s.Get(j,0), p );
    const R norm = Pow( sum, 1/p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
SchattenNorm( const DistMatrix<F,U,V>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("SchattenNorm");
#endif
    typedef typename Base<F>::type R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SVD( B, s );

    // TODO: Think of how to make this more stable
    const int kLocal = s.LocalHeight();
    R localSum = 0;
    for( int j=kLocal-1; j>=0; --j ) 
        localSum += Pow( s.GetLocal(j,0), p );
    R sum;
    mpi::AllReduce( &localSum, &sum, 1, mpi::SUM, A.Grid().VRComm() );
    const R norm = Pow( sum, 1/p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V>
inline typename Base<F>::type
HermitianSchattenNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("HermitianSchattenNorm");
#endif
    typedef typename Base<F>::type R;

    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, B, s );

    // TODO: Think of how to make this more stable
    const int kLocal = s.LocalHeight();
    R localSum = 0;
    for( int j=kLocal-1; j>=0; --j )
        localSum += Pow( s.GetLocal(j,0), p );
    R sum;
    mpi::AllReduce( &localSum, &sum, 1, mpi::SUM, A.Grid().VRComm() );
    const R norm = Pow( sum, 1/p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V>
inline typename Base<F>::type
SymmetricSchattenNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("SymmetricSchattenNorm");
#endif
    typedef typename Base<F>::type R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVD( B, s );

    // TODO: Think of how to make this more stable
    const int kLocal = s.LocalHeight();
    R localSum = 0;
    for( int j=0; j<kLocal; ++j )
        localSum += Pow( s.GetLocal(j,0), p );
    R sum;
    mpi::AllReduce( &localSum, &sum, 1, mpi::SUM, A.Grid().VRComm() );
    const R norm = Pow( sum, 1/p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_SCHATTEN_HPP
