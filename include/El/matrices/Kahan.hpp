/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_KAHAN_HPP
#define EL_KAHAN_HPP

// I haven't decided on the appropriate generalization to complex cosine/sine
// pairs. For now, given phi, we will compute the corresponding partner as the
// real value sqrt(1-|phi|^2)

namespace El {

template<typename F> 
inline void
MakeKahan( Matrix<F>& A, F phi )
{
    DEBUG_ONLY(CallStackEntry cse("MakeKahan"))
    typedef Base<F> Real;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= Real(1) || Abs(phi) == Real(0) )
        LogicError("|phi| must be in (0,1)");

    const F zeta = Sqrt(F(1)-phi*Conj(phi));

    MakeZeros( A );
    for( Int i=0; i<n; ++i )
    {
        const F zetaPow = Pow( zeta, Real(i) );
        A.Set( i, i, zetaPow );
        for( Int j=1; j<n; ++j )
            A.Set( i, j, -phi*zetaPow );
    }
}

template<typename F>
inline void
MakeKahan( AbstractDistMatrix<F>& A, F phi )
{
    DEBUG_ONLY(CallStackEntry cse("MakeKahan"))
    typedef Base<F> Real;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= Real(1) || Abs(phi) == Real(0) )
        LogicError("|phi| must be in (0,1)");

    const F zeta = Sqrt(F(1)-phi*Conj(phi));

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        const F zetaPow = Pow( zeta, Real(i) );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
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
MakeKahan( AbstractBlockDistMatrix<F>& A, F phi )
{
    DEBUG_ONLY(CallStackEntry cse("MakeKahan"))
    typedef Base<F> Real;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Kahan");
    if( Abs(phi) >= Real(1) || Abs(phi) == Real(0) )
        LogicError("|phi| must be in (0,1)");

    const F zeta = Sqrt(F(1)-phi*Conj(phi));

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        const F zetaPow = Pow( zeta, Real(i) );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
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
    DEBUG_ONLY(CallStackEntry cse("Kahan"))
    A.Resize( n, n );
    MakeKahan( A, phi );
}

template<typename F>
inline void
Kahan( AbstractDistMatrix<F>& A, Int n, F phi )
{
    DEBUG_ONLY(CallStackEntry cse("Kahan"))
    A.Resize( n, n );
    MakeKahan( A, phi );
}

template<typename F>
inline void
Kahan( AbstractBlockDistMatrix<F>& A, Int n, F phi )
{
    DEBUG_ONLY(CallStackEntry cse("Kahan"))
    A.Resize( n, n );
    MakeKahan( A, phi );
}

} // namespace El

#endif // ifndef EL_KAHAN_HPP
