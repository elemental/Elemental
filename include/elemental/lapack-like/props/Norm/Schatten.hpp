/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_NORM_SCHATTEN_HPP
#define ELEM_NORM_SCHATTEN_HPP

#include ELEM_MAKESYMMETRIC_INC
#include ELEM_MAKEHERMITIAN_INC
#include ELEM_SVD_INC

namespace elem {

template<typename F> 
inline Base<F>
SchattenNorm( const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SchattenNorm"))
    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    SVD( B, s );

    // TODO: Think of how to make this more stable
    const Int k = s.Height();
    R sum = 0;
    for( Int j=k-1; j>=0; --j )
        sum += Pow( s.Get(j,0), p ); 
    return Pow( sum, 1/p ); 
}

template<typename F>
inline Base<F>
HermitianSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSchattenNorm"))
    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    HermitianSVD( uplo, B, s );

    // TODO: Think of how to make this more stable
    const Int k = s.Height();
    R sum = 0;
    for( Int j=k-1; j>=0; --j )
        sum += Pow( s.Get(j,0), p );
    return Pow( sum, 1/p );
}

template<typename F>
inline Base<F>
SymmetricSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSchattenNorm"))
    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    MakeSymmetric( uplo, B );
    SVD( B, s );

    // TODO: Think of how to make this more stable
    const Int k = s.Height();
    R sum = 0;
    for( Int j=0; j<k; ++j )
        sum += Pow( s.Get(j,0), p );
    return Pow( sum, 1/p );
}

template<typename F,Dist U,Dist V> 
inline Base<F>
SchattenNorm( const DistMatrix<F,U,V>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SchattenNorm"))
    typedef Base<F> R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SVD( B, s );

    // TODO: Think of how to make this more stable
    const Int kLocal = s.LocalHeight();
    R localSum = 0;
    for( Int j=kLocal-1; j>=0; --j ) 
        localSum += Pow( s.GetLocal(j,0), p );
    const R sum = mpi::AllReduce( localSum, A.Grid().VRComm() );
    return Pow( sum, 1/p );
}

template<typename F,Dist U,Dist V>
inline Base<F>
HermitianSchattenNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSchattenNorm"))
    typedef Base<F> R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, B, s );

    // TODO: Think of how to make this more stable
    const Int kLocal = s.LocalHeight();
    R localSum = 0;
    for( Int j=kLocal-1; j>=0; --j )
        localSum += Pow( s.GetLocal(j,0), p );
    const R sum = mpi::AllReduce( localSum, A.Grid().VRComm() );
    return Pow( sum, 1/p );
}

template<typename F,Dist U,Dist V>
inline Base<F>
SymmetricSchattenNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSchattenNorm"))
    typedef Base<F> R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVD( B, s );

    // TODO: Think of how to make this more stable
    const Int kLocal = s.LocalHeight();
    R localSum = 0;
    for( Int j=0; j<kLocal; ++j )
        localSum += Pow( s.GetLocal(j,0), p );
    const R sum = mpi::AllReduce( localSum, A.Grid().VRComm() );
    return Pow( sum, 1/p );
}

} // namespace elem

#endif // ifndef ELEM_NORM_SCHATTEN_HPP
