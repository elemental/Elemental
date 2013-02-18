/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANNORM_TWO_HPP
#define LAPACK_HERMITIANNORM_TWO_HPP

#include "elemental/blas-like/level1/MakeHermitian.hpp"

#ifndef WITHOUT_PMRRR
#include "elemental/lapack-like/HermitianSVD.hpp"
#endif // ifndef WITHOUT_PMRRR
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/SVD.hpp"

namespace elem {

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
// TODO: Enable support for sequential MRRR
/*
#ifndef WITHOUT_PMRRR
    HermitianSingularValues( uplo, B, s );
#else
    MakeHermitian( uplo, B );
    SingularValues( B, s );
#endif
*/
    MakeHermitian( uplo, B );
    SingularValues( B, s );
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
#ifndef WITHOUT_PMRRR
    HermitianSingularValues( uplo, B, s );
#else
    MakeHermitian( uplo, B );
    SingularValues( B, s );
#endif
    const R norm = InfinityNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_HERMITIANNORM_TWO_HPP
