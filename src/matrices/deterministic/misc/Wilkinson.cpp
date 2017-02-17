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
void Wilkinson( Matrix<T>& A, Int k )
{
    EL_DEBUG_CSE
    const Int n = 2*k+1;
    Zeros( A, n, n );
    FillDiagonal( A, T(1), -1 );
    FillDiagonal( A, T(1),  1 );
    
    for( Int j=0; j<=k; ++j )
        A.Set( j, j, T(k-j) );
    for( Int j=k+1; j<n; ++j )
        A.Set( j, j, T(j-k) );
}

template<typename T>
void Wilkinson( AbstractDistMatrix<T>& A, Int k )
{
    EL_DEBUG_CSE
    const Int n = 2*k+1;
    Zeros( A, n, n );
    FillDiagonal( A, T(1), -1 );
    FillDiagonal( A, T(1),  1 );
    
    for( Int j=0; j<=k; ++j )
        A.Set( j, j, T(k-j) );
    for( Int j=k+1; j<n; ++j )
        A.Set( j, j, T(j-k) );
}

#define PROTO(T) \
  template void Wilkinson( Matrix<T>& A, Int k ); \
  template void Wilkinson( AbstractDistMatrix<T>& A, Int k );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
