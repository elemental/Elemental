/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SKELETON_HPP
#define ELEM_LAPACK_SKELETON_HPP

#include "elemental/blas-like/level1/Adjoint.hpp"
#include "elemental/lapack-like/QR/BusingerGolub.hpp"
#include "elemental/lapack-like/Pseudoinverse.hpp"

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
  Matrix<Int>& pR, Matrix<Int>& pC, 
  Matrix<F>& Z, Int maxSteps, BASE(F) tol )
{
    DEBUG_ONLY(CallStackEntry cse("Skeleton"))
    // Find the row permutation
    Matrix<F> B;
    Adjoint( A, B );
    qr::BusingerGolub( B, pR, maxSteps, tol );
    const Int numSteps = pR.Height();

    // Form pinv(AR')=pinv(AR)'
    Adjoint( A, B );
    ApplyColumnPivots( B, pR );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form K := A pinv(AR)
    Matrix<F> K;
    Gemm( NORMAL, ADJOINT, F(1), A, B, K );

    // Find the column permutation (force the same number of steps)
    B = A;
    qr::BusingerGolub( B, pC, numSteps );

    // Form pinv(AC)
    B = A;
    ApplyColumnPivots( B, pC );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form Z := pinv(AC) K = pinv(AC) (A pinv(AR))
    Gemm( NORMAL, NORMAL, F(1), B, K, Z );
}

template<typename F> 
inline void
Skeleton
( const DistMatrix<F>& A, 
  DistMatrix<Int,VR,STAR>& pR, DistMatrix<Int,VR,STAR>& pC, 
  DistMatrix<F>& Z, Int maxSteps, BASE(F) tol )
{
    DEBUG_ONLY(CallStackEntry cse("Skeleton"))
    const Grid& g = A.Grid();

    // Find the row permutation
    DistMatrix<F> B(g);
    Adjoint( A, B );
    qr::BusingerGolub( B, pR, maxSteps, tol );
    const Int numSteps = pR.Height();

    // Form pinv(AR')=pinv(AR)'
    Adjoint( A, B );
    ApplyColumnPivots( B, pR );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form K := A pinv(AR)
    DistMatrix<F> K;
    Gemm( NORMAL, ADJOINT, F(1), A, B, K );

    // Find the column permutation (force the same number of steps)
    B = A;
    qr::BusingerGolub( B, pC, numSteps );

    // Form pinv(AC)
    B = A;
    ApplyColumnPivots( B, pC );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form Z := pinv(AC) K = pinv(AC) (A pinv(AR))
    Gemm( NORMAL, NORMAL, F(1), B, K, Z );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_SKELETON_HPP
