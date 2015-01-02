/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T,typename S>
void Shift( Matrix<T>& A, S alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Shift"))
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            A.Update( i, j, alpha );
}

template<typename T,typename S>
void Shift( AbstractDistMatrix<T>& A, S alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Shift"))
    Shift( A.Matrix(), alpha );
}

template<typename T,typename S>
void Shift( AbstractBlockDistMatrix<T>& A, S alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Shift"))
    Shift( A.Matrix(), alpha );
}

template<typename T,typename S>
void Shift( DistMultiVec<T>& A, S alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Shift"))
    Shift( A.Matrix(), alpha );
}

#define PROTO_TYPES(T,S) \
  template void Shift( Matrix<T>& A, S alpha ); \
  template void Shift( AbstractDistMatrix<T>& A, S alpha ); \
  template void Shift( AbstractBlockDistMatrix<T>& A, S alpha ); \
  template void Shift( DistMultiVec<T>& A, S alpha );

#define PROTO_SAME(T) PROTO_TYPES(T,T) \

#define PROTO_INT(T) PROTO_SAME(T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_SAME(T)

#define PROTO_COMPLEX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_SAME(T)

#include "El/macros/Instantiate.h"

} // namespace El
