/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LEGENDRE_HPP
#define ELEM_LEGENDRE_HPP

#include ELEM_ZEROS_INC

namespace elem {

template<typename F> 
inline void
MakeLegendre( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeLegendre"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

    const Int n = A.Width();
    for( Int j=0; j<n-1; ++j )
    {
        const F gamma = F(1) / Pow( F(2)*F(j+1), F(2) );
        const F beta = F(1) / (F(2)*Sqrt(F(1)-gamma));
        A.Set( j+1, j, beta );
        A.Set( j, j+1, beta );
    }
}

template<typename F,Dist U,Dist V>
inline void
MakeLegendre( DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeLegendre"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            if( j == i+1 || j == i-1 )
            {
                const Int k = Max( i, j );
                const F gamma = F(1) / Pow( F(2)*F(k), F(2) );
                const F beta = F(1) / (F(2)*Sqrt(F(1)-gamma));
                A.SetLocal( iLoc, jLoc, beta );
            }
        }
    }
}

template<typename F,Dist U,Dist V>
inline void
MakeLegendre( BlockDistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeLegendre"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            if( j == i+1 || j == i-1 )
            {
                const Int k = Max( i, j );
                const F gamma = F(1) / Pow( F(2)*F(k), F(2) );
                const F beta = F(1) / (F(2)*Sqrt(F(1)-gamma));
                A.SetLocal( iLoc, jLoc, beta );
            }
        }
    }
}

template<typename F> 
inline void
Legendre( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Legendre"))
    A.Resize( n, n );
    MakeLegendre( A );
}

template<typename F,Dist U,Dist V> 
inline void
Legendre( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Legendre"))
    A.Resize( n, n );
    MakeLegendre( A );
}

template<typename F,Dist U,Dist V> 
inline void
Legendre( BlockDistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Legendre"))
    A.Resize( n, n );
    MakeLegendre( A );
}

} // namespace elem

#endif // ifndef ELEM_LEGENDRE_HPP
