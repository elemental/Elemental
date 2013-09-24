/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_KAHAN_HPP
#define ELEM_MATRICES_KAHAN_HPP

// I haven't decided on the appropriate generalization to complex cosine/sine
// pairs. For now, given phi, we will compute the corresponding partner as the
// real value sqrt(1-|phi|^2)

namespace elem {

template<typename F> 
inline void
MakeKahan( Matrix<F>& A, F phi )
{
#ifndef RELEASE
    CallStackEntry cse("MakeKahan");
#endif
    typedef BASE(F) R;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= R(1) || Abs(phi) == R(0) )
        LogicError("|phi| must be in (0,1)");

    const F zeta = Sqrt(1-phi*Conj(phi));

    MakeZeros( A );
    for( Int i=0; i<n; ++i )
    {
        const F zetaPow = Pow( zeta, R(i) );
        A.Set( i, i, zetaPow );
        for( Int j=1; j<n; ++j )
            A.Set( i, j, -phi*zetaPow );
    }
}

template<typename F,Distribution U,Distribution V>
inline void
MakeKahan( DistMatrix<F,U,V>& A, F phi )
{
#ifndef RELEASE
    CallStackEntry cse("MakeKahan");
#endif
    typedef BASE(F) R;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= R(1) || Abs(phi) == R(0) )
        LogicError("|phi| must be in (0,1)");

    const F zeta = Sqrt(1-phi*Conj(phi));

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = colShift + iLoc*colStride;
        const F zetaPow = Pow( zeta, R(i) );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            if( i > j )       
                A.SetLocal( iLoc, jLoc, F(0) ); 
            else if( i == j )
                A.SetLocal( iLoc, jLoc, zetaPow );
            else
                A.SetLocal( iLoc, jLoc, -phi*zetaPow );
        }
    }
}

template<typename F>
inline void
Kahan( Matrix<F>& A, Int n, F phi )
{
#ifndef RELEASE
    CallStackEntry cse("Kahan");
#endif
    A.ResizeTo( n, n );
    MakeKahan( A, phi );
}

template<typename F>
inline Matrix<F>
Kahan( Int n, F phi )
{
    Matrix<F> A( n, n );
    MakeKahan( A, phi );
    return A;
}

template<typename F,Distribution U,Distribution V>
inline void
Kahan( DistMatrix<F,U,V>& A, Int n, F phi )
{
#ifndef RELEASE
    CallStackEntry cse("Kahan");
#endif
    A.ResizeTo( n, n );
    MakeKahan( A, phi );
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Kahan( const Grid& g, Int n, F phi )
{
    DistMatrix<F,U,V> A( n, n, g );
    MakeKahan( A, phi );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_KAHAN_HPP
