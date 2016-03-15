/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO: Generalize to support several columns at once?
template<typename F>
bool LatticeCoordinates
( const Matrix<F>& B,
  const Matrix<F>& y,
        Matrix<F>& x ) 
{
    DEBUG_ONLY(CSE cse("LatticeCoordinates"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    if( y.Height() != m || y.Width() != 1 )
        LogicError("y should have been an ",m," x 1 vector");

    if( FrobeniusNorm(y) == Real(0) )
    {
        Zeros( x, n, 1 );
        return true;
    }
    
    Matrix<F> BRed( B );
    Matrix<F> UB, RB;
    auto infoB = LLL( BRed, UB, RB );
    auto MB = BRed( ALL, IR(0,infoB.rank) );

    Matrix<F> A;
    Zeros( A, m, infoB.rank+1 );
    {
        auto AL = A( ALL, IR(0,infoB.rank) ); 
        auto aR = A( ALL, IR(infoB.rank) );
        AL = MB;
        aR = y;
    }
    // Reduce A in-place
    Matrix<F> UA, RA;
    auto infoA = LLL( A, UA, RA );
    if( infoA.nullity != 1 )
        return false;

    // Solve for x_M such that M_B x_M = y
    // NOTE: The last column of U_A should hold the coordinates of the single
    //       member of the null-space of (the original) A
    Matrix<F> xM;
    xM = UA( IR(0,infoA.rank), IR(infoB.rank) );
    if( UA.Get(infoA.rank,infoB.rank) == F(1) )
        xM *= F(-1);
    DEBUG_ONLY(
      else if( UA.Get(infoA.rank,infoB.rank) != F(-1) )
          RuntimeError("Invalid member of null space");
    )

    // Map xM back to the original coordinates using the portion of the 
    // unimodular transformation of B (U_B) which produced the image of B 
    auto UBM = UB( ALL, IR(0,infoB.rank) );
    Zeros( x, n, 1 );
    Gemv( NORMAL, F(1), UBM, xM, F(0), x );
    
    /*
    if( infoB.nullity != 0 )
    {
        Matrix<F> C;
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

#define PROTO(F) \
  template bool LatticeCoordinates \
  ( const Matrix<F>& B, \
    const Matrix<F>& y, \
          Matrix<F>& x );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
