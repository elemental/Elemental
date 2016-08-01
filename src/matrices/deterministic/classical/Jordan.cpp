/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/matrices.hpp>

namespace El {

template<typename T>
void Jordan( Matrix<T>& J, Int n, T lambda )
{
    DEBUG_CSE
    Zeros( J, n, n );
    FillDiagonal( J, lambda );
    FillDiagonal( J, T(1), 1 );
}

template<typename T>
void Jordan( AbstractDistMatrix<T>& J, Int n, T lambda )
{
    DEBUG_CSE
    Zeros( J, n, n );
    FillDiagonal( J, lambda );
    FillDiagonal( J, T(1), 1 );
}

#define PROTO(T) \
  template void Jordan( Matrix<T>& J, Int n, T lambda ); \
  template void Jordan( AbstractDistMatrix<T>& J, Int n, T lambda );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
