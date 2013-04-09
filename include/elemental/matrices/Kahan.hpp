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
Kahan( F phi, int n, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Kahan");
#endif
    A.ResizeTo( n, n );
    MakeKahan( phi, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
Kahan( F phi, int n, DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Kahan");
#endif
    A.ResizeTo( n, n );
    MakeKahan( phi, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
MakeKahan( F phi, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("MakeKahan");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= R(1) )
        throw std::logic_error("Phi must be in (0,1)");

    const F zeta = Sqrt(1-phi*Conj(phi));

    MakeZeros( A );
    for( int i=0; i<n; ++i )
    {
        const F zetaPow = Pow( zeta, R(i) );
        A.Set( i, i, zetaPow );
        for( int j=1; j<n; ++j )
            A.Set( i, j, -phi*zetaPow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
MakeKahan( F phi, DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeKahan");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= R(1) )
        throw std::logic_error("Phi must be in (0,1)");

    const F zeta = Sqrt(1-phi*Conj(phi));

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = colShift + iLocal*colStride;
        const F zetaPow = Pow( zeta, R(i) );
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const int j = rowShift + jLocal*rowStride;
            if( i > j )       
                A.SetLocal( iLocal, jLocal, F(0) ); 
            else if( i == j )
                A.SetLocal( iLocal, jLocal, zetaPow );
            else
                A.SetLocal( iLocal, jLocal, -phi*zetaPow );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef MATRICES_KAHAN_HPP
