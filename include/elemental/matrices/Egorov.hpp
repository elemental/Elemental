/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_EGOROV_HPP
#define ELEM_MATRICES_EGOROV_HPP

namespace elem {

template<typename Real,class RealFunctor> 
inline void
MakeEgorov( Matrix<Complex<Real> >& A, const RealFunctor& phase )
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

template<typename Real,Distribution U,Distribution V,class RealFunctor>
inline void
MakeEgorov( DistMatrix<Complex<Real>,U,V>& A, const RealFunctor& phase )
{
    DEBUG_ONLY(CallStackEntry cse("MakeEgorov"))
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            const Real theta = phase(i,j);
            const Real realPart = Cos(theta);
            const Real imagPart = Sin(theta);
            A.SetLocal( iLoc, jLoc, Complex<Real>(realPart,imagPart) );
        }
    }
}

template<typename Real,class RealFunctor>
inline void
Egorov( Matrix<Complex<Real> >& A, const RealFunctor& phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    MakeEgorov( A, phase );
}

template<typename Real,class RealFunctor>
inline Matrix<Complex<Real> >
Egorov( const RealFunctor& phase, Int n )
{
    Matrix<Complex<Real>> A( n, n );
    MakeEgorov( A, phase );
    return A;
}

template<typename Real,Distribution U,Distribution V,class RealFunctor>
inline void
Egorov( DistMatrix<Complex<Real>,U,V>& A, const RealFunctor& phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    MakeEgorov( A, phase );
}

template<typename Real,Distribution U=MC,Distribution V=MR,class RealFunctor>
inline DistMatrix<Complex<Real>,U,V>
Egorov( const Grid& g, const RealFunctor& phase, Int n )
{
    DistMatrix<Complex<Real>,U,V> A( n, n, g );
    MakeEgorov( A, phase );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_EGOROV_HPP
