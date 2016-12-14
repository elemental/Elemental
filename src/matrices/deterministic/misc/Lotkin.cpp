/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/matrices.hpp>

namespace El {

template<typename F>
void Lotkin( Matrix<F>& A, Int n )
{
    EL_DEBUG_CSE
    Hilbert( A, n );
    // Set first row to all ones
    for( Int j=0; j<n; ++j )
        A.Set( 0, j, F(1) );
}

template<typename F>
void Lotkin( AbstractDistMatrix<F>& A, Int n )
{
    EL_DEBUG_CSE
    Hilbert( A, n );
    // Set first row to all ones
    if( A.ColShift() == 0 )
    {
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            A.SetLocal( 0, jLoc, F(1) );
    } 
}

#define PROTO(F) \
  template void Lotkin( Matrix<F>& A, Int n ); \
  template void Lotkin( AbstractDistMatrix<F>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
