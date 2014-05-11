/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HILBERT_HPP
#define ELEM_HILBERT_HPP

namespace elem {

template<typename F> 
inline void
MakeHilbert( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHilbert"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Hilbert");

    const F one = F(1);
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, one/F(i+j+1) );
}

template<typename F,Dist U,Dist V>
inline void
MakeHilbert( DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHilbert"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Hilbert");

    const F one = F(1);
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, one/F(i+j+1) );
        }
    }
}

template<typename F,Dist U,Dist V>
inline void
MakeHilbert( BlockDistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHilbert"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Hilbert");

    const F one = F(1);
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, one/F(i+j+1) );
        }
    }
}

template<typename F>
inline void
Hilbert( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Hilbert"))
    A.Resize( n, n );
    MakeHilbert( A );
}

template<typename F,Dist U,Dist V>
inline void
Hilbert( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Hilbert"))
    A.Resize( n, n );
    MakeHilbert( A );
}

template<typename F,Dist U,Dist V>
inline void
Hilbert( BlockDistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Hilbert"))
    A.Resize( n, n );
    MakeHilbert( A );
}

} // namespace elem

#endif // ifndef ELEM_HILBERT_HPP
