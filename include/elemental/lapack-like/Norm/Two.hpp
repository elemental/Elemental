/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_TWO_HPP
#define LAPACK_NORM_TWO_HPP

#include "elemental/blas-like/level1/MakeSymmetric.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"
#include "elemental/lapack-like/SVD.hpp"

namespace elem {

template<typename F> 
inline typename Base<F>::type
TwoNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TwoNorm");
#endif
    typedef typename Base<F>::type R;
    Matrix<F> B( A );
    Matrix<R> s;
    SVD( B, s );
    const R norm = InfinityNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianTwoNorm");
#endif
    typedef typename Base<F>::type R;
    Matrix<F> B( A );
    Matrix<R> s;
    HermitianSVD( uplo, B, s );
    const R norm = InfinityNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricTwoNorm");
#endif
    typedef typename Base<F>::type R;
    Matrix<F> B( A );
    Matrix<R> s;
    MakeSymmetric( uplo, B );
    SVD( B, s );
    const R norm = MaxNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
TwoNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("TwoNorm");
#endif
    typedef typename Base<F>::type R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SVD( B, s );
    const R norm = InfinityNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V>
inline typename Base<F>::type
HermitianTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianTwoNorm");
#endif
    typedef typename Base<F>::type R;
    DistMatrix<F,U,V> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, B, s );
    const R norm = InfinityNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V>
inline typename Base<F>::type
SymmetricTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricTwoNorm");
#endif
    typedef typename Base<F>::type R;
    DistMatrix<F,U,V> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVD( B, s );
    const R norm = MaxNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_TWO_HPP
