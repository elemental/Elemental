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
T Dot( const Matrix<T>& A, const Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dot"))
    return HilbertSchmidt( A, B );
}

template<typename T>
T Dot( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dot"))
    return HilbertSchmidt( A, B );
}

template<typename T>
T Dot( const DistMultiVec<T>& A, const DistMultiVec<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dot"))
    return HilbertSchmidt( A, B );
}

#define PROTO(T) \
  template T Dot( const Matrix<T>& A, const Matrix<T>& B ); \
  template T Dot \
  ( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B ); \
  template T Dot( const DistMultiVec<T>& A, const DistMultiVec<T>& B );

#include "El/macros/Instantiate.h"

} // namespace El
