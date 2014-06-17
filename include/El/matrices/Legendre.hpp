/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LEGENDRE_HPP
#define EL_LEGENDRE_HPP



namespace El {

template<typename F> 
inline void
MakeLegendre( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeLegendre"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Legendre");
    Zero( A );

    const Int n = A.Width();
    for( Int j=0; j<n-1; ++j )
    {
        const F gamma = F(1) / Pow( F(2)*F(j+1), F(2) );
        const F beta = F(1) / (F(2)*Sqrt(F(1)-gamma));
        A.Set( j+1, j, beta );
        A.Set( j, j+1, beta );
    }
}

template<typename F>
inline void
MakeLegendre( AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeLegendre"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Legendre");
    Zero( A );

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
MakeLegendre( AbstractBlockDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeLegendre"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Legendre");
    Zero( A );

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

template<typename F> 
inline void
Legendre( AbstractDistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Legendre"))
    A.Resize( n, n );
    MakeLegendre( A );
}

template<typename F> 
inline void
Legendre( AbstractBlockDistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Legendre"))
    A.Resize( n, n );
    MakeLegendre( A );
}

} // namespace El

#endif // ifndef EL_LEGENDRE_HPP
