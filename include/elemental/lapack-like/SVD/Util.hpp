/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SVD_UTIL_HPP
#define ELEM_LAPACK_SVD_UTIL_HPP

#include "elemental/blas-like/level1/Adjoint.hpp"
#include "elemental/lapack-like/Norm/One.hpp"

namespace elem {
namespace svd {

template<typename F>
inline void
CheckScale( DistMatrix<F>& A, bool& needRescaling, BASE(F)& scale )
{
    typedef BASE(F) R;

    scale = 1;
    needRescaling = false;
    const R oneNormOfA = OneNorm( A );
    const R safeMin = lapack::MachineSafeMin<R>();
    const R precision = lapack::MachinePrecision<R>();
    const R smallNumber = safeMin/precision;
    const R bigNumber = 1/smallNumber;
    const R rhoMin = Sqrt(smallNumber);
    const R rhoMax = Min( Sqrt(bigNumber), 1/Sqrt(Sqrt(safeMin)) );

    if( oneNormOfA > 0 && oneNormOfA < rhoMin )
    {
        needRescaling = true;
        scale = rhoMin/oneNormOfA;
    }
    else if( oneNormOfA > rhoMax )
    {
        needRescaling = true;
        scale = rhoMax/oneNormOfA;
    }
}

template<typename F>
inline void
DivideAndConquerSVD( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V )
{
#ifndef RELEASE
    CallStackEntry entry("svd::DivideAndConquerSVD");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min(m,n);
    s.ResizeTo( k, 1 );
    Matrix<F> U( m, k );
    Matrix<F> VAdj( k, n );
    lapack::DivideAndConquerSVD
    ( m, n, A.Buffer(), A.LDim(), s.Buffer(), U.Buffer(), U.LDim(),
      VAdj.Buffer(), VAdj.LDim() );

    A = U;
    Adjoint( VAdj, V );
}

template<typename F>
inline void
QRSVD( Matrix<F>& A, Matrix<BASE(F)>& s, Matrix<F>& V )
{
#ifndef RELEASE
    CallStackEntry entry("svd::QRSVD");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min(m,n);
    s.ResizeTo( k, 1 );
    Matrix<F> U( m, k );
    Matrix<F> VAdj( k, n );
    lapack::QRSVD
    ( m, n, A.Buffer(), A.LDim(), s.Buffer(), U.Buffer(), U.LDim(),
      VAdj.Buffer(), VAdj.LDim() );

    A = U;
    Adjoint( VAdj, V );
}

} // namespace svd
} // namespace elem

#endif // ifndef ELEM_LAPACK_SVD_UTIL_HPP
