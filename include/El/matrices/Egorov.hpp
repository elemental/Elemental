/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_EGOROV_HPP
#define EL_EGOROV_HPP

namespace El {

template<typename Real,class RealFunctor> 
inline void
MakeEgorov( Matrix<Complex<Real>>& A, const RealFunctor& phase )
{
    DEBUG_ONLY(CallStackEntry cse("MakeEgorov"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real theta = phase(i,j);
            const Real realPart = Cos(theta);
            const Real imagPart = Sin(theta);
            A.Set( i, j, Complex<Real>(realPart,imagPart) );
        }
    }
}

template<typename Real,class RealFunctor>
inline void
MakeEgorov( AbstractDistMatrix<Complex<Real>>& A, const RealFunctor& phase )
{
    DEBUG_ONLY(CallStackEntry cse("MakeEgorov"))
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            const Real theta = phase(i,j);
            const Real realPart = Cos(theta);
            const Real imagPart = Sin(theta);
            A.SetLocal( iLoc, jLoc, Complex<Real>(realPart,imagPart) );
        }
    }
}

template<typename Real,class RealFunctor>
inline void
MakeEgorov
( AbstractBlockDistMatrix<Complex<Real>>& A, const RealFunctor& phase )
{
    DEBUG_ONLY(CallStackEntry cse("MakeEgorov"))
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            const Real theta = phase(i,j);
            const Real realPart = Cos(theta);
            const Real imagPart = Sin(theta);
            A.SetLocal( iLoc, jLoc, Complex<Real>(realPart,imagPart) );
        }
    }
}

template<typename Real,class RealFunctor>
inline void
Egorov( Matrix<Complex<Real>>& A, const RealFunctor& phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    MakeEgorov( A, phase );
}

template<typename Real,class RealFunctor>
inline void
Egorov( AbstractDistMatrix<Complex<Real>>& A, const RealFunctor& phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    MakeEgorov( A, phase );
}

template<typename Real,class RealFunctor>
inline void
Egorov
( AbstractBlockDistMatrix<Complex<Real>>& A, const RealFunctor& phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    MakeEgorov( A, phase );
}

} // namespace El

#endif // ifndef EL_EGOROV_HPP
