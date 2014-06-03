/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IDENTITY_HPP
#define EL_IDENTITY_HPP

namespace El {

template<typename T> 
inline void
MakeIdentity( Matrix<T>& I )
{
    DEBUG_ONLY(CallStackEntry cse("MakeIdentity"))
    Zero( I );
    const Int m = I.Height();
    const Int n = I.Width();
    for( Int j=0; j<std::min(m,n); ++j )
        I.Set( j, j, T(1) );
}

template<typename T>
inline void
MakeIdentity( AbstractDistMatrix<T>& I )
{
    DEBUG_ONLY(CallStackEntry cse("MakeIdentity"))
    Zero( I.Matrix() );

    const Int localHeight = I.LocalHeight();
    const Int localWidth = I.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = I.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = I.GlobalRow(iLoc);
            if( i == j )
                I.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

template<typename T>
inline void
MakeIdentity( AbstractBlockDistMatrix<T>& I )
{
    DEBUG_ONLY(CallStackEntry cse("MakeIdentity"))
    Zero( I.Matrix() );

    const Int localHeight = I.LocalHeight();
    const Int localWidth = I.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = I.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = I.GlobalRow(iLoc);
            if( i == j )
                I.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

template<typename T>
inline void
Identity( Matrix<T>& I, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Identity"))
    I.Resize( m, n );
    MakeIdentity( I );
}

template<typename T>
inline void
Identity( AbstractDistMatrix<T>& I, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Identity"))
    I.Resize( m, n );
    MakeIdentity( I );
}

template<typename T>
inline void
Identity( AbstractBlockDistMatrix<T>& I, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Identity"))
    I.Resize( m, n );
    MakeIdentity( I );
}

} // namespace El

#endif // ifndef EL_IDENTITY_HPP
