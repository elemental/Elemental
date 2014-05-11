/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LYAPUNOV_HPP
#define ELEM_LYAPUNOV_HPP

#include ELEM_ADJOINT_INC
#include "./Sylvester.hpp"

namespace elem {

// A is assumed to have all of its eigenvalues in the open right-half plane.
// X is then returned as the solution of the system of equations:
//    A X + X A^H = C
//
// See Chapter 2 of Nicholas J. Higham's "Functions of Matrices"

template<typename F>
inline void
Lyapunov
( const Matrix<F>& A, const Matrix<F>& C, Matrix<F>& X, 
  SignCtrl<Base<F>> signCtrl=SignCtrl<Base<F>>() )
{
    DEBUG_ONLY(
        CallStackEntry cse("Lyapunov");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( C.Height() != A.Height() || C.Width() != A.Height() )
            LogicError("C must conform with A");
    )
    const Int m = A.Height();
    Matrix<F> W, WTL, WTR,
                 WBL, WBR;
    Zeros( W, 2*m, 2*m );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    WTL = A;
    Adjoint( A, WBR ); Scale( F(-1), WBR );
    WTR = C; Scale( F(-1), WTR );
    Sylvester( m, W, X, signCtrl );
}

template<typename F>
inline void
Lyapunov
( const DistMatrix<F>& A, const DistMatrix<F>& C, DistMatrix<F>& X,
  SignCtrl<Base<F>> signCtrl=SignCtrl<Base<F>>() )
{
    DEBUG_ONLY(
        CallStackEntry cse("Sylvester");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( C.Height() != A.Height() || C.Width() != A.Height() )
            LogicError("C must conform with A");
        if( A.Grid() != C.Grid() )
            LogicError("A and C must have the same grid");
    )
    const Grid& g = A.Grid();
    const Int m = A.Height();
    DistMatrix<F> W(g), WTL(g), WTR(g),
                        WBL(g), WBR(g);
    Zeros( W, 2*m, 2*m );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    WTL = A;
    Adjoint( A, WBR ); Scale( F(-1), WBR );
    WTR = C; Scale( F(-1), WTR );
    Sylvester( m, W, X, signCtrl );
}

} // namespace elem

#endif // ifndef ELEM_LYAPUNOV_HPP
