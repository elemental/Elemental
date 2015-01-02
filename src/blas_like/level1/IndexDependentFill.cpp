/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void IndexDependentFill( Matrix<T>& A, std::function<T(Int,Int)> func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentFill"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func(i,j) );
}

template<typename T>
void IndexDependentFill
( AbstractDistMatrix<T>& A, std::function<T(Int,Int)> func )
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

template<typename T>
void IndexDependentFill
( AbstractBlockDistMatrix<T>& A, std::function<T(Int,Int)> func )
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

#define PROTO(T) \
  template void IndexDependentFill \
  ( Matrix<T>& A, std::function<T(Int,Int)> func ); \
  template void IndexDependentFill \
  ( AbstractDistMatrix<T>& A, std::function<T(Int,Int)> func ); \
  template void IndexDependentFill \
  ( AbstractBlockDistMatrix<T>& A, std::function<T(Int,Int)> func );

#include "El/macros/Instantiate.h"

} // namespace El
