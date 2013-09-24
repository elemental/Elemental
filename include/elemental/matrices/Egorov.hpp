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
MakeEgorov( Matrix<Complex<R> >& A, const RealFunctor& phase )
{
#ifndef RELEASE
    CallStackEntry cse("MakeEgorov");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
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
    CallStackEntry cse("MakeEgorov");
#endif
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
            const R theta = phase(i,j);
            const R realPart = cos(theta);
            const R imagPart = sin(theta);
            A.SetLocal( iLoc, jLoc, Complex<R>(realPart,imagPart) );
        }
    }
}

template<typename R,class RealFunctor>
inline void
Egorov( Matrix<Complex<R> >& A, const RealFunctor& phase, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Egorov");
#endif
    A.ResizeTo( n, n );
    MakeEgorov( A, phase );
}

template<typename R,class RealFunctor>
inline Matrix<Complex<R> >
Egorov( const RealFunctor& phase, Int n )
{
    Matrix<Complex<R> > A( n, n );
    MakeEgorov( A, phase );
    return A;
}

template<typename R,Distribution U,Distribution V,class RealFunctor>
inline void
Egorov( DistMatrix<Complex<R>,U,V>& A, const RealFunctor& phase, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Egorov");
#endif
    A.ResizeTo( n, n );
    MakeEgorov( A, phase );
}

template<typename R,Distribution U=MC,Distribution V=MR,class RealFunctor>
inline DistMatrix<Complex<R>,U,V>
Egorov( const Grid& g, const RealFunctor& phase, Int n )
{
    DistMatrix<Complex<R>,U,V> A( n, n, g );
    MakeEgorov( A, phase );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_EGOROV_HPP
