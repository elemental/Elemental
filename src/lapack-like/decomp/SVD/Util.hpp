/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SVD_UTIL_HPP
#define EL_SVD_UTIL_HPP

namespace El {
namespace svd {

template<typename F>
inline bool
CheckScale( DistMatrix<F>& A, Base<F>& scale )
{
    scale = 1;
    typedef Base<F> Real;
    const Real oneNormOfA = OneNorm( A );
    const Real safeMin = lapack::MachineSafeMin<Real>();
    const Real precision = lapack::MachinePrecision<Real>();
    const Real smallNumber = safeMin/precision;
    const Real bigNumber = 1/smallNumber;
    const Real rhoMin = Sqrt(smallNumber);
    const Real rhoMax = Min( Sqrt(bigNumber), 1/Sqrt(Sqrt(safeMin)) );

    if( oneNormOfA > 0 && oneNormOfA < rhoMin )
    {
        scale = rhoMin/oneNormOfA;
        return true;
    }
    else if( oneNormOfA > rhoMax )
    {
        scale = rhoMax/oneNormOfA;
        return true;
    }
    else
        return false;
}

template<typename F>
inline void
DivideAndConquerSVD( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("svd::DivideAndConquerSVD"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min(m,n);
    s.Resize( k, 1 );
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
QRSVD( Matrix<F>& A, Matrix<Base<F>>& s, Matrix<F>& V )
{
    DEBUG_ONLY(CallStackEntry cse("svd::QRSVD"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min(m,n);
    s.Resize( k, 1 );
    Matrix<F> U( m, k );
    Matrix<F> VAdj( k, n );
    lapack::QRSVD
    ( m, n, A.Buffer(), A.LDim(), s.Buffer(), U.Buffer(), U.LDim(),
      VAdj.Buffer(), VAdj.LDim() );

    A = U;
    Adjoint( VAdj, V );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_UTIL_HPP
