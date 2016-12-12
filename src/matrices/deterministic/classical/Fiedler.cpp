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

template<typename Field>
void Fiedler( Matrix<Field>& A, const vector<Field>& c )
{
    EL_DEBUG_CSE
    const Int n = c.size();
    A.Resize( n, n );
    auto fiedlerFill = [&]( Int i, Int j ) { return Abs(c[i]-c[j]); };
    IndexDependentFill( A, function<Field(Int,Int)>(fiedlerFill) );
}

template<typename Field>
void Fiedler( AbstractDistMatrix<Field>& A, const vector<Field>& c )
{
    EL_DEBUG_CSE
    const Int n = c.size();
    A.Resize( n, n );
    auto fiedlerFill = [&]( Int i, Int j ) { return Abs(c[i]-c[j]); };
    IndexDependentFill( A, function<Field(Int,Int)>(fiedlerFill) );
}

#define PROTO(Field) \
  template void Fiedler( Matrix<Field>& A, const vector<Field>& c ); \
  template void Fiedler( AbstractDistMatrix<Field>& A, const vector<Field>& c );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
