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
void EntrywiseFill( Matrix<T>& A, std::function<T(void)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseFill"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func() );
}

template<typename T>
void EntrywiseFill( AbstractDistMatrix<T>& A, std::function<T(void)> func )
{ EntrywiseFill( A.Matrix(), func ); }

template<typename T>
void EntrywiseFill( AbstractBlockDistMatrix<T>& A, std::function<T(void)> func )
{ EntrywiseFill( A.Matrix(), func ); }

template<typename T>
void EntrywiseFill( DistMultiVec<T>& A, std::function<T(void)> func )
{ EntrywiseFill( A.Matrix(), func ); }

#define PROTO(T) \
  template void EntrywiseFill \
  ( Matrix<T>& A, std::function<T(void)> func ); \
  template void EntrywiseFill \
  ( AbstractDistMatrix<T>& A, std::function<T(void)> func ); \
  template void EntrywiseFill \
  ( AbstractBlockDistMatrix<T>& A, std::function<T(void)> func ); \
  template void EntrywiseFill \
  ( DistMultiVec<T>& A, std::function<T(void)> func );

#include "El/macros/Instantiate.h"

} // namespace El
