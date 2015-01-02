/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace axpy {
namespace util {

template<typename T>
void InterleaveMatrixUpdate
( T alpha, Int height, Int width,
  const T* A, Int colStrideA, Int rowStrideA,
        T* B, Int colStrideB, Int rowStrideB )
{
    // TODO: Add OpenMP parallelization and/or optimize
    for( Int j=0; j<width; ++j )
        blas::Axpy
        ( height, alpha,
          &A[rowStrideA*j], colStrideA,
          &B[rowStrideB*j], colStrideB );
}

#define PROTO(T) \
  template void InterleaveMatrixUpdate \
  ( T alpha, Int height, Int width, \
    const T* A, Int colStrideA, Int rowStrideA, \
          T* B, Int colStrideB, Int rowStrideB );

#include "El/macros/Instantiate.h"

} // namespace util
} // namespace axpy
} // namespace El
