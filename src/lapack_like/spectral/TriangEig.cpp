/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2015-2016, Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "./TriangEig/MultiShiftSolve.hpp"

namespace El {

template<typename Field>
void TriangEig( Matrix<Field>& U, Matrix<Field>& X )
{

    EL_DEBUG_CSE
    const Int m = U.Height();

    // Make X the negative of the strictly upper triangle of  U
    X = U;
    MakeTrapezoidal( UPPER, X, 1 );
    Scale( Field(-1), X );

    // Solve multi-shift triangular system
    Matrix<Field> shifts, scales;
    GetDiagonal( U, shifts );
    // The following is a specialized alternative to
    //  SafeMultiShiftTrsm
    //  ( LEFT, UPPER, NORMAL, Field(1), U, shifts, X, scales );
    triang_eig::MultiShiftSolve( U, shifts, X, scales );
    SetDiagonal( X, scales );

    // Normalize eigenvectors
    for( Int j=0; j<m; ++j )
    {
        auto xj = X( IR(0,j+1), IR(j) );
        Scale( 1/Nrm2(xj), xj );
    }
}

template<typename Field>
void TriangEig
( const AbstractDistMatrix<Field>& UPre,
        AbstractDistMatrix<Field>& XPre )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Field,Field,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<Field,Field,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    // Make X the negative of the strictly upper triangle of  U
    X = U;
    MakeTrapezoidal( UPPER, X, 1 );
    Scale( Field(-1), X );

    // Solve multi-shift triangular system
    const Grid& g = U.Grid();
    DistMatrix<Field,VR,STAR> shifts(g), scales(g);
    GetDiagonal( U, shifts );
    // The following is a specialized alternative to
    //  SafeMultiShiftTrsm
    //  ( LEFT, UPPER, NORMAL, Field(1), U, shifts, X, scales );
    triang_eig::MultiShiftSolve( U, shifts, X, scales );
    SetDiagonal( X, scales );

    // Normalize eigenvectors
    // TODO(poulson): Exploit the upper-triangular structure
    DistMatrix<Base<Field>,MR,STAR> colNorms(g);
    ColumnTwoNorms( X, colNorms );
    DiagonalSolve( RIGHT, NORMAL, colNorms, X );
}

#define PROTO(Field) \
  template void TriangEig \
  (       Matrix<Field>& T, \
          Matrix<Field>& X ); \
  template void TriangEig \
  ( const AbstractDistMatrix<Field>& T, \
          AbstractDistMatrix<Field>& X );

#define EL_NO_INT_PROTO
#define EL_ENABLE_QUAD
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
