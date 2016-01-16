/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// NOTE: There are *many* algorithms for (pseudo-)skeleton/CUR decompositions,
//       and, for now, we will simply implement one.

// TODO: More algorithms and more options (e.g., default tolerances).

// TODO: Implement randomized algorithms from Jiawei Chiu and Laurent Demanet's 
//       "Sublinear randomized algorithms for skeleton decompositions"?

namespace El {

template<typename F> 
void Skeleton
( const Matrix<F>& A, 
        Permutation& PR,
        Permutation& PC,
        Matrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("Skeleton"))
    // Find the row permutation
    Matrix<F> B;
    Adjoint( A, B );
    Matrix<F> t;
    Matrix<Base<F>> d;
    QR( B, t, d, PR, ctrl );
    const Int numSteps = t.Height();

    // Form pinv(AR')=pinv(AR)'
    Adjoint( A, B );
    PR.PermuteCols( B );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form K := A pinv(AR)
    Matrix<F> K;
    Gemm( NORMAL, ADJOINT, F(1), A, B, K );

    // Find the column permutation (force the same number of steps)
    B = A;
    auto secondCtrl = ctrl; 
    secondCtrl.adaptive = false;
    secondCtrl.boundRank = true;
    secondCtrl.maxRank = numSteps;
    QR( B, t, d, PC, secondCtrl );

    // Form pinv(AC)
    B = A;
    PC.PermuteCols( B );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form Z := pinv(AC) K = pinv(AC) (A pinv(AR))
    Gemm( NORMAL, NORMAL, F(1), B, K, Z );
}

template<typename F> 
void Skeleton
( const ElementalMatrix<F>& APre, 
        DistPermutation& PR,
        DistPermutation& PC,
        ElementalMatrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("Skeleton"))

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    const Grid& g = A.Grid();

    // Find the row permutation
    DistMatrix<F> B(g);
    Adjoint( A, B );
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    QR( B, t, d, PR, ctrl );
    const Int numSteps = t.Height();

    // Form pinv(AR')=pinv(AR)'
    Adjoint( A, B );
    PR.PermuteCols( B );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form K := A pinv(AR)
    DistMatrix<F> K(g);
    Gemm( NORMAL, ADJOINT, F(1), A, B, K );

    // Find the column permutation (force the same number of steps)
    B = A;
    auto secondCtrl = ctrl; 
    secondCtrl.adaptive = false;
    secondCtrl.boundRank = true;
    secondCtrl.maxRank = numSteps;
    QR( B, t, d, PC, secondCtrl );

    // Form pinv(AC)
    B = A;
    PC.PermuteCols( B );
    B.Resize( B.Height(), numSteps );
    Pseudoinverse( B );

    // Form Z := pinv(AC) K = pinv(AC) (A pinv(AR))
    Gemm( NORMAL, NORMAL, F(1), B, K, Z );
}

#define PROTO(F) \
  template void Skeleton \
  ( const Matrix<F>& A, \
          Permutation& PR, \
          Permutation& PC, \
          Matrix<F>& Z, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void Skeleton \
  ( const ElementalMatrix<F>& A, \
          DistPermutation& PR, \
          DistPermutation& PC, \
          ElementalMatrix<F>& Z, \
    const QRCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
// NOTE: These cannot be enabled until there is more general SVD support
/*
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
*/
#include "El/macros/Instantiate.h"

} // namespace El
