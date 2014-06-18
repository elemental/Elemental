/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_INDEXDEPENDENTFILL_HPP
#define EL_INDEXDEPENDENTFILL_HPP

namespace El {

template<typename T,class Function>
inline void
IndexDependentFill( Matrix<T>& A, Function func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentFill"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func(i,j) );
}

template<typename T,class Function>
inline void
IndexDependentFill( AbstractDistMatrix<T>& A, Function func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentFill"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, func(i,j) );
        }
    }
}

template<typename T,class Function>
inline void
IndexDependentFill( AbstractBlockDistMatrix<T>& A, Function func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentFill"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, func(i,j) );
        }
    }
}

} // namespace El

#endif // ifndef EL_INDEXDEPENDENTFILL_HPP
