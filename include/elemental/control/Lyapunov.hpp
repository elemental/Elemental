/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONTROL_LYAPUNOV_HPP
#define ELEM_CONTROL_LYAPUNOV_HPP

#include "elemental/blas-like/level1/Adjoint.hpp"
#include "elemental/control/Sylvester.hpp"

namespace elem {

// A is assumed to have all of its eigenvalues in the open right-half plane.
// X is then returned as the solution of the system of equations:
//    A X + X A^H = C
//
// See Chapter 2 of Nicholas J. Higham's "Functions of Matrices"

template<typename F>
inline int
Lyapunov( const Matrix<F>& A, const Matrix<F>& C, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Lyapunov");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( C.Height() != A.Height() || C.Width() != A.Height() )
        LogicError("C must conform with A");
#endif
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
    return Sylvester( m, W, X );
}

template<typename F>
inline int
Lyapunov( const DistMatrix<F>& A, const DistMatrix<F>& C, DistMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Sylvester");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( C.Height() != A.Height() || C.Width() != A.Height() )
        LogicError("C must conform with A");
    if( A.Grid() != C.Grid() )
        LogicError("A and C must have the same grid");
#endif
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
    return Sylvester( m, W, X );
}

} // namespace elem

#endif // ifndef ELEM_CONTROL_LYAPUNOV_HPP
