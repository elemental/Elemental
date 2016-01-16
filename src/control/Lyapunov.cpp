/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// A is assumed to have all of its eigenvalues in the open right-half plane.
// X is then returned as the solution of the system of equations:
//    A X + X A^H = C
//
// See Chapter 2 of Nicholas J. Higham's "Functions of Matrices"

template<typename F>
void Lyapunov
( const Matrix<F>& A, const Matrix<F>& C, Matrix<F>& X, 
  SignCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(
      CSE cse("Lyapunov");
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
    Adjoint( A, WBR ); WBR *= -1;
    WTR = C;           WTR *= -1;
    Sylvester( m, W, X, ctrl );
}

template<typename F>
void Lyapunov
( const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& C, 
        ElementalMatrix<F>& X, SignCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(
      CSE cse("Sylvester");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( C.Height() != A.Height() || C.Width() != A.Height() )
          LogicError("C must conform with A");
      AssertSameGrids( A, C );
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
    Adjoint( A, WBR ); WBR *= -1;
    WTR = C;           WTR *= -1;
    Sylvester( m, W, X, ctrl );
}

#define PROTO(F) \
  template void Lyapunov \
  ( const Matrix<F>& A, \
    const Matrix<F>& C, \
          Matrix<F>& X, \
    SignCtrl<Base<F>> ctrl ); \
  template void Lyapunov \
  ( const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& C, \
          ElementalMatrix<F>& X, \
    SignCtrl<Base<F>> ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
