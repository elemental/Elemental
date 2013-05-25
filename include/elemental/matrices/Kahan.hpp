/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_KAHAN_HPP
#define MATRICES_KAHAN_HPP

// I haven't decided on the appropriate generalization to complex cosine/sine
// pairs. For now, given phi, we will compute the corresponding partner as the
// real value sqrt(1-|phi|^2)

namespace elem {

template<typename F>
inline void
Kahan( Matrix<F>& A, int n, F phi )
{
#ifndef RELEASE
    CallStackEntry entry("Kahan");
#endif
    A.ResizeTo( n, n );
    MakeKahan( A, phi );
}

template<typename F,Distribution U,Distribution V>
inline void
Kahan( DistMatrix<F,U,V>& A, int n, F phi )
{
#ifndef RELEASE
    CallStackEntry entry("Kahan");
#endif
    A.ResizeTo( n, n );
    MakeKahan( A, phi );
}

template<typename F> 
inline void
MakeKahan( Matrix<F>& A, F phi )
{
#ifndef RELEASE
    CallStackEntry entry("MakeKahan");
#endif
    typedef BASE(F) R;

    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= R(1) || Abs(phi) == R(0) )
        throw std::logic_error("|phi| must be in (0,1)");

    const F zeta = Sqrt(1-phi*Conj(phi));

    MakeZeros( A );
    for( int i=0; i<n; ++i )
    {
        const F zetaPow = Pow( zeta, R(i) );
        A.Set( i, i, zetaPow );
        for( int j=1; j<n; ++j )
            A.Set( i, j, -phi*zetaPow );
    }
}

template<typename F,Distribution U,Distribution V>
inline void
MakeKahan( DistMatrix<F,U,V>& A, F phi )
{
#ifndef RELEASE
    CallStackEntry entry("MakeKahan");
#endif
    typedef BASE(F) R;

    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= R(1) || Abs(phi) == R(0) )
        throw std::logic_error("|phi| must be in (0,1)");

    const F zeta = Sqrt(1-phi*Conj(phi));

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*colStride;
        const F zetaPow = Pow( zeta, R(i) );
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*rowStride;
            if( i > j )       
                A.SetLocal( iLoc, jLoc, F(0) ); 
            else if( i == j )
                A.SetLocal( iLoc, jLoc, zetaPow );
            else
                A.SetLocal( iLoc, jLoc, -phi*zetaPow );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_KAHAN_HPP
