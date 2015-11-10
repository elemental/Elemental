/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO: Sparse matrix versions

template<typename T>
void Round( Matrix<T>& A )
{
    DEBUG_ONLY(CSE cse("Round"))
    const Int m = A.Height();
    const Int n = A.Width();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            ABuf[i+j*ALDim] = Round(ABuf[i+j*ALDim]);
}

template<>
void Round( Matrix<Int>& A )
{ }

template<typename T>
void Round( AbstractDistMatrix<T>& A )
{ Round( A.Matrix() ); }

template<typename T>
void Round( DistMultiVec<T>& A )
{ Round( A.Matrix() ); }

#define PROTO(T) \
  template void Round( Matrix<T>& A ); \
  template void Round( AbstractDistMatrix<T>& A ); \
  template void Round( DistMultiVec<T>& A );

#include "El/macros/Instantiate.h"

} // namespace El
