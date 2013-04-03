/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP
#define CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/lapack-like/SVD.hpp"
#include "elemental/convex/SoftThreshold.hpp"

namespace elem {

template<typename F>
inline int
SingularValueSoftThreshold( Matrix<F>& A, typename Base<F>::type tau )
{
#ifndef RELEASE
    PushCallStack("SingularValueSoftThreshold");
#endif
    typedef typename Base<F>::type R;
    Matrix<F> U( A );
    Matrix<R> s;
    Matrix<F> V;

    // TODO: Exploit zeros in soft-thresholded singular values (with custom SVD)
    SVD( U, s, V );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    const int rank = ZeroNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return rank;
}

template<typename F>
inline int
SingularValueSoftThreshold( DistMatrix<F>& A, typename Base<F>::type tau )
{
#ifndef RELEASE
    PushCallStack("SingularValueSoftThreshold");
#endif
    typedef typename Base<F>::type R;
    DistMatrix<F> U( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    DistMatrix<F> V( A.Grid() );

    svd::Thresholded( U, s, V, tau );
    SoftThreshold( s, tau );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    const int rank = ZeroNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return rank;
}

} // namespace elem

#endif // ifndef CONVEX_SINGULARVALUESOFTTHRESHOLD_HPP
