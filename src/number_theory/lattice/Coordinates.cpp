/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2016, Ron Estrin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// TODO(poulson): Generalize to support several columns at once?
template<typename Field>
bool LatticeCoordinates
( const Matrix<Field>& B,
  const Matrix<Field>& y,
        Matrix<Field>& x )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    if( y.Height() != m || y.Width() != 1 )
        LogicError("y should have been an ",m," x 1 vector");

    if( FrobeniusNorm(y) == Real(0) )
    {
        Zeros( x, n, 1 );
        return true;
    }

    Matrix<Field> BRed( B );
    Matrix<Field> UB, RB;
    auto infoB = LLL( BRed, UB, RB );
    auto MB = BRed( ALL, IR(0,infoB.rank) );

    Matrix<Field> A;
    Zeros( A, m, infoB.rank+1 );
    {
        auto AL = A( ALL, IR(0,infoB.rank) );
        auto aR = A( ALL, IR(infoB.rank) );
        AL = MB;
        aR = y;
    }
    // Reduce A in-place
    Matrix<Field> UA, RA;
    auto infoA = LLL( A, UA, RA );
    if( infoA.nullity != 1 )
        return false;

    // Solve for x_M such that M_B x_M = y
    // NOTE: The last column of U_A should hold the coordinates of the single
    //       member of the null-space of (the original) A
    Matrix<Field> xM;
    xM = UA( IR(0,infoA.rank), IR(infoB.rank) );
    const Field gamma = UA(infoA.rank,infoB.rank);
    if( Abs(gamma) != Real(1) )
        LogicError("Invalid member of null space");
    else
        xM *= -Conj(gamma);

    // Map xM back to the original coordinates using the portion of the
    // unimodular transformation of B (U_B) which produced the image of B
    auto UBM = UB( ALL, IR(0,infoB.rank) );
    Zeros( x, n, 1 );
    Gemv( NORMAL, Field(1), UBM, xM, Field(0), x );

    /*
    if( infoB.nullity != 0 )
    {
        Matrix<Field> C;
        Zeros( C, m, infoB.nullity+1 );
        auto cL = C( ALL, IR(infoB.rank-1) );
        auto CR = C( ALL, IR(infoB.rank,END) );

        // Reduce the kernel of B
        CR = UB( ALL, IR(infoB.rank,END) );
        LLL( CR );

        // Attempt to reduce the (reduced) kernel out of the coordinates
        // TODO: Which column to grab from the result?!?
        cL = x;
        LLL( C );
        x = cL;
    }
    */

    return true;
}

#define PROTO(Field) \
  template bool LatticeCoordinates \
  ( const Matrix<Field>& B, \
    const Matrix<Field>& y, \
          Matrix<Field>& x );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
