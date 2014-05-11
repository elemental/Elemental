/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QR_SOLVEAFTER_HPP
#define ELEM_QR_SOLVEAFTER_HPP

#include ELEM_ZERO_INC
#include ELEM_TRSM_INC
#include ELEM_QR_INC

// TODO: Extend for BusingerGolub support

namespace elem {
namespace qr {

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
    DEBUG_ONLY(CallStackEntry cse("qr::SolveAfter"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m < n )
        LogicError("Must have full column rank");
    // TODO: Add scaling
    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            LogicError("A and B do not conform");

        // Copy B into X
        X = B;

        // Apply Q' to X
        qr::ApplyQ( LEFT, ADJOINT, A, t, d, X );

        // Shrink X to its new height
        X.Resize( n, X.Width() );

        // Solve against R (checking for singularities)
        auto AT = LockedView( A, 0, 0, n, n );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AT, X, true );
    }
    else // orientation in {TRANSPOSE,ADJOINT}
    {
        if( n != B.Height() )
            LogicError("A and B do not conform");

        // Copy B into X
        X.Resize( m, B.Width() );
        Matrix<F> XT, XB;
        PartitionDown( X, XT, XB, n );
        XT = B;
        Zero( XB );

        if( orientation == TRANSPOSE )
            Conjugate( XT );

        // Solve against R' (checking for singularities)
        auto AT = LockedView( A, 0, 0, n, n );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AT, XT, true );

        // Apply Q to X
        qr::ApplyQ( LEFT, NORMAL, A, t, d, X );

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
    DEBUG_ONLY(CallStackEntry cse("qr::SolveAfter"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    if( m < n )
        LogicError("Must have full column rank");

    // TODO: Add scaling

    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            LogicError("A and B do not conform");

        // Copy B into X
        X = B;

        if( orientation == TRANSPOSE )
            Conjugate( X );

        // Apply Q' to X
        qr::ApplyQ( LEFT, ADJOINT, A, t, d, X );

        // Shrink X to its new height
        X.Resize( n, X.Width() );

        // Solve against R (checking for singularities)
        auto AT = LockedView( A, 0, 0, n, n );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AT, X, true );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
    else
    {
        // Copy B into X
        X.Resize( m, B.Width() );
        DistMatrix<F> XT(g), XB(g);
        PartitionDown( X, XT, XB, n );
        XT = B;
        Zero( XB );

        if( orientation == TRANSPOSE )
            Conjugate( XT );

        // Solve against R' (checking for singularities)
        auto AT = LockedView( A, 0, 0, n, n );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AT, XT, true );

        // Apply Q to X
        qr::ApplyQ( LEFT, NORMAL, A, t, d, X );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_QR_SOLVEAFTER_HPP
