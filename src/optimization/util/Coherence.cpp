/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
Base<Field> Coherence( const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    Matrix<Field> B( A );
    Matrix<Base<Field>> norms;
    ColumnTwoNorms( B, norms );

    DiagonalSolve( RIGHT, NORMAL, norms, B, true );
    Matrix<Field> C;
    Identity( C, A.Width(), A.Width() );
    Herk( UPPER, ADJOINT, Base<Field>(-1), B, Base<Field>(1), C );

    return HermitianMaxNorm( UPPER, C );
}

template<typename Field>
Base<Field> Coherence( const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    DistMatrix<Field> B( A );
    DistMatrix<Base<Field>,MR,STAR> norms(B.Grid());
    ColumnTwoNorms( B, norms );

    DiagonalSolve( RIGHT, NORMAL, norms, B, true );
    DistMatrix<Field> C(B.Grid());
    Identity( C, A.Width(), A.Width() );
    Herk( UPPER, ADJOINT, Base<Field>(-1), B, Base<Field>(1), C );

    return HermitianMaxNorm( UPPER, C );
}

#define PROTO(Field) \
  template Base<Field> Coherence( const Matrix<Field>& A ); \
  template Base<Field> Coherence( const AbstractDistMatrix<Field>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
