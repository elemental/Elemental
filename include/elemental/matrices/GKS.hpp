/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GKS_HPP
#define ELEM_GKS_HPP

// The Golub Klema Stewart matrix is upper-triangular with 1/sqrt(j) on its 
// j'th diagonal entry and -1/sqrt(j) elsewhere in the upper triangle.
// 
// It was originally introduced as an example of where greedy RRQR fails.

namespace elem {

template<typename F> 
inline void
MakeGKS( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGKS"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix GKS");

    MakeZeros( A );
    for( Int j=0; j<n; ++j )
    {
        const F jDiag = F(1)/Sqrt(F(j+1));
        for( Int i=0; i<j; ++i )
            A.Set( i, j, -jDiag );
        A.Set( j, j, jDiag );
    }
}

template<typename F,Dist U,Dist V>
inline void
MakeGKS( DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGKS"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix GKS");

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const F jDiag = F(1)/Sqrt(F(j+1));
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            if( i < j )
                A.SetLocal( iLoc, jLoc, -jDiag );
            else if( i == j )
                A.SetLocal( iLoc, jLoc, jDiag );
            else
                A.SetLocal( iLoc, jLoc, 0 );
        }
    }
}

template<typename F,Dist U,Dist V>
inline void
MakeGKS( BlockDistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGKS"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix GKS");

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const F jDiag = F(1)/Sqrt(F(j+1));
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            if( i < j )
                A.SetLocal( iLoc, jLoc, -jDiag );
            else if( i == j )
                A.SetLocal( iLoc, jLoc, jDiag );
            else
                A.SetLocal( iLoc, jLoc, 0 );
        }
    }
}

template<typename F>
inline void
GKS( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GKS"))
    A.Resize( n, n );
    MakeGKS( A );
}

template<typename F,Dist U,Dist V>
inline void
GKS( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GKS"))
    A.Resize( n, n );
    MakeGKS( A );
}

template<typename F,Dist U,Dist V>
inline void
GKS( BlockDistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GKS"))
    A.Resize( n, n );
    MakeGKS( A );
}

} // namespace elem

#endif // ifndef ELEM_GKS_HPP
