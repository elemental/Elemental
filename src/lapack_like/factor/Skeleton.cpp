/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// TODO: Use column-pivoted LQ decompositions to avoid unnecessary explicit
//       adjoints

// TODO: Avoid reapplying Q when forming K

// NOTE: There are *many* algorithms for (pseudo-)skeleton/CUR decompositions,
//       and, for now, we will simply implement one.

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
    DEBUG_CSE
    Matrix<F> AAdj;
    Adjoint( A, AAdj );

    // Find the row permutation
    Matrix<F> B(AAdj);
    Matrix<F> phase;
    Matrix<Base<F>> signature;
    QR( B, phase, signature, PR, ctrl );
    const Int numSteps = phase.Height();
    B.Resize( B.Height(), numSteps );
    // Form K' := (A pinv(AR))' = pinv(AR') A'
    Matrix<F> KAdj;
    qr::SolveAfter( NORMAL, B, phase, signature, AAdj, KAdj );
    // Form K := (K')'
    Matrix<F> K;
    Adjoint( KAdj, K );

    // Find the column permutation (force the same number of steps)
    B = A;
    auto secondCtrl = ctrl; 
    secondCtrl.adaptive = false;
    secondCtrl.boundRank = true;
    secondCtrl.maxRank = numSteps;
    QR( B, phase, signature, PC, secondCtrl );
    // Form Z := pinv(AC) K = pinv(AC) (A pinv(AR))
    B.Resize( B.Height(), numSteps );
    qr::SolveAfter( NORMAL, B, phase, signature, K, Z );
}

template<typename F> 
void Skeleton
( const ElementalMatrix<F>& APre, 
        DistPermutation& PR,
        DistPermutation& PC,
        ElementalMatrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();
    const Grid& g = A.Grid();

    DistMatrix<F> AAdj(g);
    Adjoint( A, AAdj );

    // Find the row permutation
    DistMatrix<F> B(AAdj);
    DistMatrix<F,MD,STAR> phase(g);
    DistMatrix<Base<F>,MD,STAR> signature(g);
    QR( B, phase, signature, PR, ctrl );
    const Int numSteps = phase.Height();
    B.Resize( B.Height(), numSteps );
    // Form K' := (A pinv(AR))' = pinv(AR') A'
    DistMatrix<F> KAdj(g);
    qr::SolveAfter( NORMAL, B, phase, signature, AAdj, KAdj );
    // Form K := (K')'
    DistMatrix<F> K(g);
    Adjoint( KAdj, K );

    // Find the column permutation (force the same number of steps)
    B = A;
    auto secondCtrl = ctrl; 
    secondCtrl.adaptive = false;
    secondCtrl.boundRank = true;
    secondCtrl.maxRank = numSteps;
    QR( B, phase, signature, PC, secondCtrl );
    // Form Z := pinv(AC) K = pinv(AC) (A pinv(AR))
    B.Resize( B.Height(), numSteps );
    qr::SolveAfter( NORMAL, B, phase, signature, K, Z );
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
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
