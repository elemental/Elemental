/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SKELETON_HPP
#define ELEM_SKELETON_HPP

#include ELEM_ADJOINT_INC
#include ELEM_QR_BUSINGERGOLUB_INC
#include ELEM_PERMUTECOLS_INC
#include ELEM_PSEUDOINVERSE_INC

// NOTE: There are *many* algorithms for (pseudo-)skeleton/CUR decompositions,
//       and, for now, we will simply implement one.

// TODO: More algorithms and more options (e.g., default tolerances).

// TODO: Implement randomized algorithms from Jiawei Chiu and Laurent Demanet's 
//       "Sublinear randomized algorithms for skeleton decompositions"?

namespace elem {

template<typename F> 
inline void
Skeleton
( const Matrix<F>& A, 
  Matrix<Int>& permR, Matrix<Int>& permC, 
  Matrix<F>& Z, Int maxSteps, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("Skeleton"))
    // Find the row permutation
    Matrix<F> B;
    Adjoint( A, B );
    const Int numSteps = qr::BusingerGolub( B, permR, maxSteps, tol );

    // Form pinv(AR')=pinv(AR)'
    Adjoint( A, B );
    InversePermuteCols( B, permR );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form K := A pinv(AR)
    Matrix<F> K;
    Gemm( NORMAL, ADJOINT, F(1), A, B, K );

    // Find the column permutation (force the same number of steps)
    B = A;
    qr::BusingerGolub( B, permC, numSteps );

    // Form pinv(AC)
    B = A;
    InversePermuteCols( B, permC );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form Z := pinv(AC) K = pinv(AC) (A pinv(AR))
    Gemm( NORMAL, NORMAL, F(1), B, K, Z );
}

template<typename F,Dist UPerm> 
inline void
Skeleton
( const DistMatrix<F>& A, 
  DistMatrix<Int,UPerm,STAR>& permR, DistMatrix<Int,UPerm,STAR>& permC, 
  DistMatrix<F>& Z, Int maxSteps, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("Skeleton"))
    const Grid& g = A.Grid();

    // Find the row permutation
    DistMatrix<F> B(g);
    Adjoint( A, B );
    const Int numSteps = qr::BusingerGolub( B, permR, maxSteps, tol );

    // Form pinv(AR')=pinv(AR)'
    Adjoint( A, B );
    InversePermuteCols( B, permR );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form K := A pinv(AR)
    DistMatrix<F> K(g);
    Gemm( NORMAL, ADJOINT, F(1), A, B, K );

    // Find the column permutation (force the same number of steps)
    B = A;
    qr::BusingerGolub( B, permC, numSteps );

    // Form pinv(AC)
    B = A;
    InversePermuteCols( B, permC );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form Z := pinv(AC) K = pinv(AC) (A pinv(AR))
    Gemm( NORMAL, NORMAL, F(1), B, K, Z );
}

} // namespace elem

#endif // ifndef ELEM_SKELETON_HPP
