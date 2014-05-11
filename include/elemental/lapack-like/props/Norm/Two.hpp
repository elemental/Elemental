/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_NORM_TWO_HPP
#define ELEM_NORM_TWO_HPP

#include ELEM_MAKESYMMETRIC_INC

#include ELEM_INFINITYNORM_INC
#include ELEM_MAXNORM_INC

#include ELEM_SVD_INC

namespace elem {

template<typename F> 
inline Base<F>
TwoNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("TwoNorm"))
    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    SVD( B, s );
    return InfinityNorm( s );
}

template<typename F>
inline Base<F>
HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTwoNorm"))
    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    HermitianSVD( uplo, B, s );
    return InfinityNorm( s );
}

template<typename F>
inline Base<F>
SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricTwoNorm"))
    typedef Base<F> R;
    Matrix<F> B( A );
    Matrix<R> s;
    MakeSymmetric( uplo, B );
    SVD( B, s );
    return MaxNorm( s );
}

template<typename F,Dist U,Dist V> 
inline Base<F>
TwoNorm( const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("TwoNorm"))
    typedef Base<F> R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SVD( B, s );
    return InfinityNorm( s );
}

template<typename F,Dist U,Dist V>
inline Base<F>
HermitianTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTwoNorm"))
    typedef Base<F> R;
    DistMatrix<F,U,V> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, B, s );
    return InfinityNorm( s );
}

template<typename F,Dist U,Dist V>
inline Base<F>
SymmetricTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricTwoNorm"))
    typedef Base<F> R;
    DistMatrix<F,U,V> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVD( B, s );
    return MaxNorm( s );
}

} // namespace elem

#endif // ifndef ELEM_NORM_TWO_HPP
