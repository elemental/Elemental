/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LQ_SOLVEAFTER_HPP
#define ELEM_LQ_SOLVEAFTER_HPP

#include ELEM_ZERO_INC
#include ELEM_TRSM_INC
#include ELEM_LQ_INC

// TODO: Extend for BusingerGolub support

namespace elem {
namespace lq {

template<typename F> 
inline void
SolveAfter
( Orientation orientation, 
  const Matrix<F>& A, 
  const Matrix<F>& t, 
  const Matrix<Base<F>>& d, 
  const Matrix<F>& B,       
        Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("lq::SolveAfter"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m > n )
        LogicError("Must have full row rank");
    // TODO: Add scaling
    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            LogicError("A and B do not conform");

        // Copy B into X
        X.Resize( n, B.Width() );
        Matrix<F> XT, XB;
        PartitionDown( X, XT, XB, m );
        XT = B;
        Zero( XB );

        // Solve against L (checking for singularities)
        auto AL = LockedView( A, 0, 0, m, m );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), AL, XT, true );

        // Apply Q' to X 
        lq::ApplyQ( LEFT, ADJOINT, A, t, d, X );
    }
    else // orientation in {TRANSPOSE,ADJOINT}
    {
        if( n != B.Height() )
            LogicError("A and B do not conform");

        // Copy B into X
        X = B;

        if( orientation == TRANSPOSE )
            Conjugate( X );

        // Apply Q to X
        lq::ApplyQ( LEFT, NORMAL, A, t, d, X );

        // Shrink X to its new height
        X.Resize( m, X.Width() );

        // Solve against L' (check for singularities)
        auto AL = LockedView( A, 0, 0, m, m );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), AL, X, true );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
}

template<typename F>
inline void
SolveAfter
( Orientation orientation,
  const DistMatrix<F>& A, 
  const DistMatrix<F,MD,STAR>& t, 
  const DistMatrix<Base<F>,MD,STAR>& d,
  const DistMatrix<F>& B, 
        DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("lq::SolveAfter"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    if( m > n )
        LogicError("Must have full row rank");

    // TODO: Add scaling

    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            LogicError("A and B do not conform");

        // Copy B into X
        X.Resize( n, B.Width() );
        DistMatrix<F> XT(g), XB(g);
        PartitionDown( X, XT, XB, m );
        XT = B;
        Zero( XB );

        if( orientation == TRANSPOSE )
            Conjugate( XT );

        // Solve against L (checking for singularities)
        auto AL = LockedView( A, 0, 0, m, m );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), AL, XT, true );

        // Apply Q' to X 
        lq::ApplyQ( LEFT, ADJOINT, A, t, d, X );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
    else
    {
        // Copy B into X
        X = B;

        if( orientation == TRANSPOSE )
            Conjugate( X );

        // Apply Q to X
        lq::ApplyQ( LEFT, NORMAL, A, t, d, X );

        // Shrink X to its new height
        X.Resize( m, X.Width() );

        // Solve against L' (check for singularities)
        auto AL = LockedView( A, 0, 0, m, m );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), AL, X, true );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
}

} // namespace lq
} // namespace elem

#endif // ifndef ELEM_LQ_SOLVEAFTER_HPP
