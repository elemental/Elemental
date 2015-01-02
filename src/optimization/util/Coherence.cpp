/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
Base<F> Coherence( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Coherence"))
    Matrix<F> B( A );
    Matrix<Base<F>> norms;
    ColumnNorms( B, norms );

    DiagonalSolve( RIGHT, NORMAL, norms, B, true );
    Matrix<F> C;
    Identity( C, A.Width(), A.Width() );
    Herk( UPPER, ADJOINT, Base<F>(-1), B, Base<F>(1), C );

    return HermitianMaxNorm( UPPER, C );
}

template<typename F>
Base<F> Coherence( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Coherence"))
    DistMatrix<F> B( A );
    DistMatrix<Base<F>,MR,STAR> norms(B.Grid());
    ColumnNorms( B, norms );

    DiagonalSolve( RIGHT, NORMAL, norms, B, true );
    DistMatrix<F> C(B.Grid());
    Identity( C, A.Width(), A.Width() );
    Herk( UPPER, ADJOINT, Base<F>(-1), B, Base<F>(1), C );

    return HermitianMaxNorm( UPPER, C );
}

#define PROTO(F) \
  template Base<F> Coherence( const Matrix<F>& A ); \
  template Base<F> Coherence( const AbstractDistMatrix<F>& A );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
