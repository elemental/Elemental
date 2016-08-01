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
void GEPPGrowth( Matrix<T>& A, Int n )
{
    DEBUG_CSE
    Identity( A, n, n );
    if( n <= 1 )
        return;

    // Set the last column to all ones
    auto aLast = A( IR(0,n), IR(n-1,n) );
    Fill( aLast, T(1) );

    // Set the subdiagonals to -1
    for( Int j=1; j<n; ++j )
        FillDiagonal( A, T(-1), -j );
}

template<typename T>
void GEPPGrowth( ElementalMatrix<T>& A, Int n )
{
    DEBUG_CSE
    Identity( A, n, n );
    if( n <= 1 )
        return;

    // Set the last column to all ones
    unique_ptr<ElementalMatrix<T>> aLast( A.Construct(A.Grid(),A.Root()) );
    View( *aLast, A, IR(0,n), IR(n-1,n) );
    Fill( *aLast, T(1) );

    // Set the subdiagonals to -1
    for( Int j=1; j<n; ++j )
        FillDiagonal( A, T(-1), -j );
}

#define PROTO(T) \
  template void GEPPGrowth( Matrix<T>& A, Int n ); \
  template void GEPPGrowth( ElementalMatrix<T>& A, Int n );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
