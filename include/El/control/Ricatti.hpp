/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RICATTI_HPP
#define EL_RICATTI_HPP

#include EL_ZEROS_INC

namespace El {

// W = | A^H  L |, where K and L are Hermitian.
//     | K   -A |
//
// The solution, X, to the equation
//   X K X - A^H X - X A = L
// is returned, as well as the number of Newton iterations for computing sgn(W).
//
// See Chapter 2 of Nicholas J. Higham's "Functions of Matrices"

template<typename F>
inline void
Ricatti
( Matrix<F>& W, Matrix<F>& X, SignCtrl<Base<F>> signCtrl=SignCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("Ricatti"))
    Sign( W, signCtrl );
    const Int n = W.Height()/2;
    Matrix<F> WTL, WTR,
              WBL, WBR;
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    // (ML, MR) = sgn(W) - I
    UpdateDiagonal( W, F(-1) );

    // Solve for X in ML X = -MR
    Matrix<F> ML, MR;
    PartitionRight( W, ML, MR, n );
    Scale( F(-1), MR );
    LeastSquares( NORMAL, ML, MR, X );
}

template<typename F>
inline void
Ricatti
( DistMatrix<F>& W, DistMatrix<F>& X, 
  SignCtrl<Base<F>> signCtrl=SignCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("Ricatti"))
    const Grid& g = W.Grid();
    Sign( W, signCtrl );
    const Int n = W.Height()/2;
    DistMatrix<F> WTL(g), WTR(g),
                  WBL(g), WBR(g);
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    // (ML, MR) = sgn(W) - I
    UpdateDiagonal( W, F(-1) );

    // Solve for X in ML X = -MR
    DistMatrix<F> ML(g), MR(g);
    PartitionRight( W, ML, MR, n );
    Scale( F(-1), MR );
    LeastSquares( NORMAL, ML, MR, X );
}

template<typename F>
inline void
Ricatti
( UpperOrLower uplo, 
  const Matrix<F>& A, const Matrix<F>& K, const Matrix<F>& L, Matrix<F>& X,
  SignCtrl<Base<F>> signCtrl=SignCtrl<Base<F>>() )
{
    DEBUG_ONLY(
        CallStackEntry cse("Sylvester");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( K.Height() != K.Width() )
            LogicError("K must be square");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
        if( A.Height() != K.Height() || A.Height() != L.Height() )
            LogicError("A, K, and L must be the same size");
    )
    const Int n = A.Height();
    Matrix<F> W, WTL, WTR,
                 WBL, WBR;
    Zeros( W, 2*n, 2*n );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    Adjoint( A, WTL );
    WBR = A; Scale( F(-1), WBR );
    WBL = K; MakeHermitian( uplo, WBL );
    WTR = L; MakeHermitian( uplo, WTR );

    Ricatti( W, X, signCtrl );
}

template<typename F>
inline void
Ricatti
( UpperOrLower uplo, 
  const DistMatrix<F>& A, const DistMatrix<F>& K, const DistMatrix<F>& L, 
  DistMatrix<F>& X, SignCtrl<Base<F>> signCtrl=SignCtrl<Base<F>>() )
{
    DEBUG_ONLY(
        CallStackEntry cse("Sylvester");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( K.Height() != K.Width() )
            LogicError("K must be square");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
        if( A.Height() != K.Height() || A.Height() != L.Height() )
            LogicError("A, K, and L must be the same size");
        if( A.Grid() != K.Grid() || K.Grid() != L.Grid() )
            LogicError("A, K, and L must have the same grid");
    )
    const Grid& g = A.Grid();
    const Int n = A.Height();
    DistMatrix<F> W(g), WTL(g), WTR(g),
                        WBL(g), WBR(g);
    Zeros( W, 2*n, 2*n );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    Adjoint( A, WTL );
    WBR = A; Scale( F(-1), WBR );
    WBL = K; MakeHermitian( uplo, WBL );
    WTR = L; MakeHermitian( uplo, WTR );

    Ricatti( W, X, signCtrl );
}

} // namespace El

#endif // ifndef EL_RICATTI_HPP
