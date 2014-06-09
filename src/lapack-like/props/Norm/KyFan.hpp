/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_KYFAN_HPP
#define EL_NORM_KYFAN_HPP

namespace El {

template<typename F> 
Base<F> KyFanNorm( const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("KyFanNorm"))
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFan norm");

    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    SVD( B, s );

    R norm = 0;
    for( Int j=k-1; j>=0; --j )
        norm += s.Get(j,0);
    return norm;
}

template<typename F>
Base<F> HermitianKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianKyFanNorm"))
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFan norm");

    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    HermitianSVD( uplo, B, s );

    R norm = 0;
    for( Int j=k-1; j>=0; --j )
        norm += s.Get(j,0);
    return norm;
}

template<typename F>
Base<F> SymmetricKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricKyFanNorm"))
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFan norm");

    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    MakeSymmetric( uplo, B );
    SVD( B, s );

    R norm = 0;
    for( Int j=k-1; j>=0; --j )
        norm += s.Get(j,0);
    return norm;
}

template<typename F,Dist U,Dist V> 
Base<F> KyFanNorm( const DistMatrix<F,U,V>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("KyFanNorm"))
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFan norm");

    typedef Base<F> R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SVD( B, s );

    R localNorm = 0;
    DistMatrix<R,VR,STAR> sTop( A.Grid() );
    LockedView( sTop, s, 0, 0, k, 1 );
    const Int localHeight = sTop.LocalHeight();
    for( Int j=localHeight-1; j>=0; --j )
        localNorm += sTop.GetLocal(j,0);
    return mpi::AllReduce( localNorm, A.Grid().VRComm() );
}

template<typename F,Dist U,Dist V>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianKyFanNorm"))
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFan norm");

    typedef Base<F> R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, B, s );

    R localNorm = 0;
    DistMatrix<R,VR,STAR> sTop( A.Grid() );
    LockedView( sTop, s, 0, 0, k, 1 );
    const Int localHeight = sTop.LocalHeight();
    for( Int j=localHeight-1; j>=0; --j )
        localNorm += sTop.GetLocal(j,0);
    return mpi::AllReduce( localNorm, A.Grid().VRComm() );
}

template<typename F,Dist U,Dist V>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricKyFanNorm"))
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFan norm");

    typedef Base<F> R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVD( B, s );

    R localNorm = 0;
    DistMatrix<R,VR,STAR> sTop( A.Grid() );
    LockedView( sTop, s, 0, 0, k, 1 );
    const Int localHeight = sTop.LocalHeight();
    for( Int j=localHeight-1; j>=0; --j )
        localNorm += sTop.GetLocal(j,0);
    return mpi::AllReduce( localNorm, A.Grid().VRComm() );
}

} // namespace El

#endif // ifndef EL_NORM_KYFAN_HPP
