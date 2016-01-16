/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// W = | A -C |, where A is m x m, B is n x n, and both are assumed to have 
//     | 0 -B |  all of their eigenvalues in the open right-half plane.
//
// The solution, X, to the equation
//   A X + X B = C
// is returned, as well as the number of Newton iterations for computing sgn(W).
//
// See Chapter 2 of Nicholas J. Higham's "Functions of Matrices"

template<typename F>
void Sylvester
( Int m,
  Matrix<F>& W,
  Matrix<F>& X,
  SignCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CSE cse("Sylvester"))
    Sign( W, ctrl );
    Matrix<F> WTL, WTR,
              WBL, WBR;
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    // WTL and WBR should be the positive and negative identity, WBL should be 
    // zero, and WTR should be -2 X
    X = WTR;
    X *= -F(1)/F(2);

    // TODO: Think of how to probe for checks on other quadrants.
    /*
    typedef Base<F> Real; 
    UpdateDiagonal( WTL, F(-1) );
    const Real errorWTL = FrobeniusNorm( WTL );
    const Int n = W.Height() - m;
    UpdateDiagonal( WBR, F(1) );
    const Real errorWBR = FrobeniusNorm( WBR );
    const Real errorWBL = FrobeniusNorm( WBL );
    */
}

template<typename F>
void Sylvester
( Int m,
  ElementalMatrix<F>& WPre,
  ElementalMatrix<F>& X, 
  SignCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CSE cse("Sylvester"))

    DistMatrixReadProxy<F,F,MC,MR> WProx( WPre );
    auto& W = WProx.Get();

    const Grid& g = W.Grid();
    Sign( W, ctrl );
    DistMatrix<F> WTL(g), WTR(g),
                  WBL(g), WBR(g);
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    // WTL and WBR should be the positive and negative identity, WBL should be 
    // zero, and WTR should be -2 X
    Copy( WTR, X );
    X *= -F(1)/F(2);

    // TODO: Think of how to probe for checks on other quadrants.
    //       Add UpdateDiagonal routine to avoid explicit identity Axpy?
    /*
    typedef Base<F> Real; 
    UpdateDiagonal( WTL, F(-1) );
    const Real errorWTL = FrobeniusNorm( WTL );
    const Int n = W.Height() - m;
    UpdateDiagonal( WBR, F(1) );
    const Real errorWBR = FrobeniusNorm( WBR );
    const Real errorWBL = FrobeniusNorm( WBL );
    */
}

template<typename F>
void Sylvester
( const Matrix<F>& A,
  const Matrix<F>& B,
  const Matrix<F>& C,
        Matrix<F>& X,
  SignCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(
      CSE cse("Sylvester");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( B.Height() != B.Width() )
          LogicError("B must be square");
      if( C.Height() != A.Height() || C.Width() != B.Height() )
          LogicError("C must conform with A and B");
    )
    const Int m = C.Height();
    const Int n = C.Width();
    Matrix<F> W, WTL, WTR,
                 WBL, WBR;
    Zeros( W, m+n, m+n );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    WTL = A;
    WBR = B; WBR *= -1;
    WTR = C; WTR *= -1;
    Sylvester( m, W, X, ctrl );
}

template<typename F>
void Sylvester
( const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& B, 
  const ElementalMatrix<F>& C,
        ElementalMatrix<F>& X, 
  SignCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(
      CSE cse("Sylvester");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( B.Height() != B.Width() )
          LogicError("B must be square");
      if( C.Height() != A.Height() || C.Width() != B.Height() )
          LogicError("C must conform with A and B");
      AssertSameGrids( A, B, C );
    )
    const Int m = C.Height();
    const Int n = C.Width();
    const Grid& g = A.Grid();
    DistMatrix<F> W(g), WTL(g), WTR(g),
                        WBL(g), WBR(g);
    Zeros( W, m+n, m+n );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    WTL = A;
    WBR = B; WBR *= -1;
    WTR = C; WTR *= -1;
    Sylvester( m, W, X, ctrl );
}

#define PROTO(F) \
  template void Sylvester \
  ( Int m, \
    Matrix<F>& W, \
    Matrix<F>& X, \
    SignCtrl<Base<F>> ctrl ); \
  template void Sylvester \
  ( Int m, \
    ElementalMatrix<F>& W, \
    ElementalMatrix<F>& X, \
    SignCtrl<Base<F>> ctrl ); \
  template void Sylvester \
  ( const Matrix<F>& A, \
    const Matrix<F>& B, \
    const Matrix<F>& C, \
          Matrix<F>& X, \
    SignCtrl<Base<F>> ctrl ); \
  template void Sylvester \
  ( const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& B, \
    const ElementalMatrix<F>& C, \
          ElementalMatrix<F>& X, \
    SignCtrl<Base<F>> ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
