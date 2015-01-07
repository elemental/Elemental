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
void FillDiagonal( Matrix<T>& A, T alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("FillDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
    {
        const Int i = j-offset;
        if( i >= 0 && i < height )
            A.Set( i, j, alpha );
    }
}

template<typename T>
void FillDiagonal( AbstractDistMatrix<T>& A, T alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("FillDiagonal"))
    const Int height = A.Height();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Int i = j-offset;
        if( i >= 0 && i < height )
            A.Set( i, j, alpha );
    }
}

template<typename T>
void FillDiagonal( AbstractBlockDistMatrix<T>& A, T alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("FillDiagonal"))
    const Int height = A.Height();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Int i = j-offset;
        if( i >= 0 && i < height )
            A.Set( i, j, alpha );
    }
}

#define PROTO(T) \
  template void FillDiagonal( Matrix<T>& A, T alpha, Int offset ); \
  template void FillDiagonal \
  ( AbstractDistMatrix<T>& A, T alpha, Int offset ); \
  template void FillDiagonal \
  ( AbstractBlockDistMatrix<T>& A, T alpha, Int offset );

#include "El/macros/Instantiate.h"

} // namespace El
