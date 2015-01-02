/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_QR_SOLVEAFTER_HPP
#define EL_QR_SOLVEAFTER_HPP

// TODO: Extend for BusingerGolub support

namespace El {
namespace qr {

template<typename F> 
void SolveAfter
( Orientation orientation, const Matrix<F>& A, 
  const Matrix<F>& t, const Matrix<Base<F>>& d, 
  const Matrix<F>& B,       Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("qr::SolveAfter"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m < n )
        LogicError("Must have full column rank");

    // TODO: Add scaling
    auto AT = A( IR(0,n), IR(0,n) );
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
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AT, XT, true );

        // Apply Q to X
        qr::ApplyQ( LEFT, NORMAL, A, t, d, X );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
}

template<typename F>
void SolveAfter
( Orientation orientation, 
  const AbstractDistMatrix<F      >& APre, const AbstractDistMatrix<F>& t, 
  const AbstractDistMatrix<Base<F>>& d,    const AbstractDistMatrix<F>& B, 
        AbstractDistMatrix<F      >& XPre )
{
    DEBUG_ONLY(CallStackEntry cse("qr::SolveAfter"))
    const Int m = APre.Height();
    const Int n = APre.Width();
    const Grid& g = APre.Grid();
    if( m < n )
        LogicError("Must have full column rank");

    auto APtr = ReadProxy<F,MC,MR>( &APre );  auto& A = *APtr;
    auto XPtr = WriteProxy<F,MC,MR>( &XPre ); auto& X = *XPtr;

    XPre.Resize( m, B.Width() );
    // TODO: Add scaling

    auto AT = A( IR(0,n), IR(0,n) );
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
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AT, X, true );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
    else
    {
        // Copy B into X
        DistMatrix<F> XT(g), XB(g);
        PartitionDown( X, XT, XB, n );
        XT = B;
        Zero( XB );

        if( orientation == TRANSPOSE )
            Conjugate( XT );

        // Solve against R' (checking for singularities)
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AT, XT, true );

        // Apply Q to X
        qr::ApplyQ( LEFT, NORMAL, A, t, d, X );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_SOLVEAFTER_HPP
