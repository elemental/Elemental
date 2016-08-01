/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_RQ_SOLVEAFTER_HPP
#define EL_RQ_SOLVEAFTER_HPP

// TODO: Extend for BusingerGolub support

namespace El {
namespace rq {

template<typename F> 
void SolveAfter
( Orientation orientation, 
  const Matrix<F>& A,
  const Matrix<F>& phase, 
  const Matrix<Base<F>>& signature,
  const Matrix<F>& B,       
        Matrix<F>& X )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    if( m > n )
        LogicError("Must have full row rank");
    // TODO: Add scaling
    auto AR = A( IR(0,m), IR(n-m,n) );
    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            LogicError("A and B do not conform");

        // Copy B into X
        X.Resize( n, B.Width() );
        auto XT = X( IR(0,m), ALL );
        auto XB = X( IR(m,n), ALL );
        XT = B;
        Zero( XB );

        // Solve against R (checking for singularities)
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AR, XT, true );

        // Apply Q' to X 
        rq::ApplyQ( LEFT, ADJOINT, A, phase, signature, X );
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
        rq::ApplyQ( LEFT, NORMAL, A, phase, signature, X );

        // Shrink X to its new height
        X.Resize( m, X.Width() );

        // Solve against R' (check for singularities)
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AR, X, true );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
}

template<typename F>
void SolveAfter
( Orientation orientation,
  const ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& phase, 
  const ElementalMatrix<Base<F>>& signature,
  const ElementalMatrix<F>& B, 
        ElementalMatrix<F>& XPre )
{
    DEBUG_CSE
    const Int m = APre.Height();
    const Int n = APre.Width();
    if( m > n )
        LogicError("Must have full row rank");

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& A = AProx.GetLocked();
    auto& X = XProx.Get();

    X.Resize( n, B.Width() );
    // TODO: Add scaling

    auto AR = A( IR(0,m), IR(n-m,n) );
    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            LogicError("A and B do not conform");

        // Copy B into X
        auto XT = X( IR(0,m), ALL );
        auto XB = X( IR(m,n), ALL );
        XT = B;
        Zero( XB );

        if( orientation == TRANSPOSE )
            Conjugate( XT );

        // Solve against R (checking for singularities)
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AR, XT, true );

        // Apply Q' to X 
        rq::ApplyQ( LEFT, ADJOINT, A, phase, signature, X );

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
        rq::ApplyQ( LEFT, NORMAL, A, phase, signature, X );

        // Shrink X to its new height
        X.Resize( m, X.Width() );

        // Solve against R' (check for singularities)
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AR, X, true );

        if( orientation == TRANSPOSE )
            Conjugate( X );
    }
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_SOLVEAFTER_HPP
