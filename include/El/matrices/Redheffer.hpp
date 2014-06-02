/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_REDHEFFER_HPP
#define EL_REDHEFFER_HPP

namespace El {

template<typename T> 
inline void
Redheffer( Matrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Redheffer"))
    R.Resize( n, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
        {
            if( j==0 || ((j+1)%(i+1))==0 )
                R.Set( i, j, T(1) );
            else
                R.Set( i, j, T(0) );
        }
    }
}

template<typename T>
inline void
Redheffer( AbstractDistMatrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Redheffer"))
    R.Resize( n, n );
    const Int localHeight = R.LocalHeight();
    const Int localWidth = R.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = R.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = R.GlobalRow(iLoc);
            if( j==0 || ((j+1)%(i+1))==0 )
                R.SetLocal( iLoc, jLoc, T(1) );
            else
                R.SetLocal( iLoc, jLoc, T(0) );
        }
    }
}

template<typename T>
inline void
Redheffer( AbstractBlockDistMatrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Redheffer"))
    R.Resize( n, n );
    const Int localHeight = R.LocalHeight();
    const Int localWidth = R.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = R.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = R.GlobalRow(iLoc);
            if( j==0 || ((j+1)%(i+1))==0 )
                R.SetLocal( iLoc, jLoc, T(1) );
            else
                R.SetLocal( iLoc, jLoc, T(0) );
        }
    }
}

} // namespace El

#endif // ifndef EL_REDHEFFER_HPP
