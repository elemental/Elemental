/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_EGOROV_HPP
#define ELEM_MATRICES_EGOROV_HPP

namespace elem {

template<typename R,class RealFunctor>
inline void
Egorov( Matrix<Complex<R> >& A, const RealFunctor& phase, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Egorov");
#endif
    A.ResizeTo( n, n );
    MakeEgorov( A, phase );
}

template<typename R,Distribution U,Distribution V,class RealFunctor>
inline void
Egorov( DistMatrix<Complex<R>,U,V>& A, const RealFunctor& phase, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Egorov");
#endif
    A.ResizeTo( n, n );
    MakeEgorov( A, phase );
}

template<typename R,class RealFunctor> 
inline void
MakeEgorov( Matrix<Complex<R> >& A, const RealFunctor& phase )
{
#ifndef RELEASE
    CallStackEntry entry("MakeEgorov");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
    {
        for( int i=0; i<m; ++i )
        {
            const R theta = phase(i,j);
            const R realPart = cos(theta);
            const R imagPart = sin(theta);
            A.Set( i, j, Complex<R>(realPart,imagPart) );
        }
    }
}

template<typename R,Distribution U,Distribution V,class RealFunctor>
inline void
MakeEgorov( DistMatrix<Complex<R>,U,V>& A, const RealFunctor& phase )
{
#ifndef RELEASE
    CallStackEntry entry("MakeEgorov");
#endif
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            const R theta = phase(i,j);
            const R realPart = cos(theta);
            const R imagPart = sin(theta);
            A.SetLocal( iLoc, jLoc, Complex<R>(realPart,imagPart) );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_EGOROV_HPP
