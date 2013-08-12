/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONTROL_RICATTI_HPP
#define ELEM_CONTROL_RICATTI_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/lapack-like/LeastSquares.hpp"
#include "elemental/lapack-like/Sign.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

// W = | A^H  L |, where K and L are Hermitian.
//     | K   -A |
//
// The solution, X, to the equation
//   X K X - A^H X - X A = L
// is returned, as well as the number of Newton iterations for computing sgn(W).
//
// See Chapter 2 of Nicholas J. Higham's "Functions of Matrices"

template<typename F>
inline int
Ricatti( Matrix<F>& W, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Ricatti");
#endif
    const Int numIts = sign::Newton( W );
    const Int n = W.Height()/2;
    Matrix<F> WTL, WTR,
              WBL, WBR;
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    // (ML, MR) = sgn(W) - I
    // TODO: Implement and replace with 'UpdateDiagonal'
    Matrix<F> I;
    Identity( I, 2*n, 2*n );
    Axpy( F(-1), I, W );

    // Solve for X in ML X = -MR
    Matrix<F> ML, MR;
    PartitionRight( W, ML, MR, n );
    Scale( F(-1), MR );
    LeastSquares( NORMAL, ML, MR, X );

    return numIts;
}

template<typename F>
inline int
Ricatti( DistMatrix<F>& W, DistMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Ricatti");
#endif
    const Grid& g = W.Grid();
    const Int numIts = sign::Newton( W );
    const Int n = W.Height()/2;
    DistMatrix<F> WTL(g), WTR(g),
                  WBL(g), WBR(g);
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    // (ML, MR) = sgn(W) - I
    // TODO: Implement and replace with 'UpdateDiagonal'
    DistMatrix<F> I(g);
    Identity( I, 2*n, 2*n );
    Axpy( F(-1), I, W );

    // Solve for X in ML X = -MR
    DistMatrix<F> ML(g), MR(g);
    PartitionRight( W, ML, MR, n );
    Scale( F(-1), MR );
    LeastSquares( NORMAL, ML, MR, X );

    return numIts;
}

template<typename F>
inline int
Ricatti
( UpperOrLower uplo, 
  const Matrix<F>& A, const Matrix<F>& K, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Sylvester");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( K.Height() != K.Width() )
        LogicError("K must be square");
    if( L.Height() != L.Width() )
        LogicError("L must be square");
    if( A.Height() != K.Height() || A.Height() != L.Height() )
        LogicError("A, K, and L must be the same size");
#endif
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

    return Ricatti( W, X );
}

template<typename F>
inline int
Ricatti
( UpperOrLower uplo, 
  const DistMatrix<F>& A, const DistMatrix<F>& K, const DistMatrix<F>& L, 
  DistMatrix<F>& X )
{
#ifndef RELEASE
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
#endif
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

    return Ricatti( W, X );
}

} // namespace elem

#endif // ifndef ELEM_CONTROL_RICATTI_HPP
