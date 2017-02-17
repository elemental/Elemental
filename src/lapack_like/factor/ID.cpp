/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// TODO: Add detailed references to Tygert et al.'s ID package and the papers
//       "Randomized algorithms for the low-rank approximation of matrices",
//       "A randomized algorithm for principal component analysis", and
//       "On the compression of low-rank matrices"

namespace El {

namespace id {

// On output, the matrix Z contains the non-trivial portion of the interpolation
// matrix, and p contains the pivots used during the iterations of
// pivoted QR. The input matrix A is unchanged.

template<typename F>
inline void
BusingerGolub
( Matrix<F>& A,
  Permutation& Omega,
  Matrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    auto ctrlCopy = ctrl;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real eps = limits::Epsilon<Real>();
    // Demand that we will be able to apply inv(R_L) to R_R by ensuring that
    // the minimum singular value is sufficiently (relatively) large
    ctrlCopy.adaptive = true;
    if( ctrl.boundRank )
    {
        ctrlCopy.tol = Max(ctrl.tol,eps*ctrl.maxRank);
    }
    else
    {
        ctrlCopy.tol = Max(ctrl.tol,eps*Min(m,n));
    }

    // Perform the pivoted QR factorization
    Matrix<F> householderScalars;
    Matrix<Base<F>> signature;
    QR( A, householderScalars, signature, Omega, ctrlCopy );
    const Int numSteps = householderScalars.Height();

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    auto RL = A( IR(0,numSteps), IR(0,numSteps) );
    auto RR = A( IR(0,numSteps), IR(numSteps,n) );
    Z = RR;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), RL, Z );
}

template<typename F>
inline void
BusingerGolub
( AbstractDistMatrix<F>& APre,
  DistPermutation& Omega,
  AbstractDistMatrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    auto ctrlCopy = ctrl;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real eps = limits::Epsilon<Real>();
    // Demand that we will be able to apply inv(R_L) to R_R by ensuring that
    // the minimum singular value is sufficiently (relatively) large
    ctrlCopy.adaptive = true;
    if( ctrl.boundRank )
    {
        ctrlCopy.tol = Max(ctrl.tol,eps*ctrl.maxRank);
    }
    else
    {
        ctrlCopy.tol = Max(ctrl.tol,eps*Min(m,n));
    }

    // Perform an adaptive pivoted QR factorization
    DistMatrix<F,MD,STAR> householderScalars(A.Grid());
    DistMatrix<Base<F>,MD,STAR> signature(A.Grid());
    QR( A, householderScalars, signature, Omega, ctrlCopy );
    const Int numSteps = householderScalars.Height();

    auto RL = A( IR(0,numSteps), IR(0,numSteps) );
    auto RR = A( IR(0,numSteps), IR(numSteps,n) );
    Copy( RR, Z );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), RL, Z );
}

} // namespace id

template<typename F>
void ID
( const Matrix<F>& A,
        Permutation& Omega,
        Matrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    Matrix<F> B( A );
    id::BusingerGolub( B, Omega, Z, ctrl );
}

template<typename F>
void ID
(       Matrix<F>& A,
        Permutation& Omega,
        Matrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl,
        bool canOverwrite )
{
    EL_DEBUG_CSE
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    id::BusingerGolub( B, Omega, Z, ctrl );
}

template<typename F>
void ID
( const AbstractDistMatrix<F>& A,
        DistPermutation& Omega,
        AbstractDistMatrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrix<F> B( A );
    id::BusingerGolub( B, Omega, Z, ctrl );
}

template<typename F>
void ID
(       AbstractDistMatrix<F>& A,
        DistPermutation& Omega,
        AbstractDistMatrix<F>& Z,
  const QRCtrl<Base<F>>& ctrl,
        bool canOverwrite )
{
    EL_DEBUG_CSE
    if( canOverwrite )
    {
        id::BusingerGolub( A, Omega, Z, ctrl );
    }
    else
    {
        DistMatrix<F> B( A );
        id::BusingerGolub( B, Omega, Z, ctrl );
    }
}

#define PROTO(F) \
  template void ID \
  ( const Matrix<F>& A, \
          Permutation& Omega, \
          Matrix<F>& Z, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void ID \
  ( const AbstractDistMatrix<F>& A, \
          DistPermutation& Omega, \
          AbstractDistMatrix<F>& Z, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void ID \
  ( Matrix<F>& A, \
    Permutation& Omega, \
    Matrix<F>& Z, \
    const QRCtrl<Base<F>>& ctrl, \
    bool canOverwrite ); \
  template void ID \
  ( AbstractDistMatrix<F>& A, \
    DistPermutation& Omega, \
    AbstractDistMatrix<F>& Z, \
    const QRCtrl<Base<F>>& ctrl, \
    bool canOverwrite );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
