/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVD_UTIL_HPP
#define EL_SVD_UTIL_HPP

namespace El {
namespace svd {

template<typename F>
bool CheckScale( Matrix<F>& A, Base<F>& scale )
{
    scale = 1;
    typedef Base<F> Real;
    const Real oneNormOfA = OneNorm( A );
    const Real safeMin = limits::SafeMin<Real>();
    const Real precision = limits::Precision<Real>();
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
bool CheckScale( AbstractDistMatrix<F>& A, Base<F>& scale )
{
    scale = 1;
    typedef Base<F> Real;
    const Real oneNormOfA = OneNorm( A );
    const Real safeMin = limits::SafeMin<Real>();
    const Real precision = limits::Precision<Real>();
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

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_UTIL_HPP
