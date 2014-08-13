/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// TODO: Add detailed references to Tygert et al.'s ID package and the papers
//       "Randomized algorithms for the low-rank approximation of matrices", 
//       "A randomized algorithm for principal component analysis", and 
//       "On the compression of low-rank matrices"

namespace El {
namespace id {

template<typename F>
inline void
PseudoTrsm( const Matrix<F>& RL, Matrix<F>& RR, Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("id::PseudoTrsm"))
    typedef Base<F> Real;
    const Int m = RR.Height();
    const Int n = RR.Width();

    // Compute the spectral radius of the triangular matrix
    Real maxAbsEig = 0;
    for( Int i=0; i<m; ++i )
        maxAbsEig = std::max( maxAbsEig, Abs(RL.Get(i,i)) );

    for( Int i=m-1; i>=0; --i )
    {
        // Apply the pseudo-inverse of the i'th diagonal value of RL 
        const F rho = RL.Get(i,i);
        const Real rhoAbs = Abs(rho);
        if( rhoAbs >= tol*maxAbsEig )
        {
            for( Int j=0; j<n; ++j ) 
            {
                const F zeta = RR.Get(i,j);
                RR.Set(i,j,zeta/rho);
            }
        }
        else
        {
            for( Int j=0; j<n; ++j )
                RR.Set(i,j,0);
        }

        // Now update RR using an outer-product of the column of RL above the 
        // i'th diagonal with the i'th row of RR
        blas::Geru
        ( i, n, 
          F(-1), RL.LockedBuffer(0,i), 1, RR.LockedBuffer(i,0), RR.LDim(), 
                 RR.Buffer(0,0), RR.LDim() );
    }
}

// For now, assume that RL is sufficiently small and give each process a full
// copy so that we may independently apply its pseudoinverse to each column of
// RR
template<typename F>
inline void
PseudoTrsm
( const AbstractDistMatrix<F>& RLPre, AbstractDistMatrix<F>& RRPre,
  Base<F> tol )
{
    DEBUG_ONLY(CallStackEntry cse("id::PseudoTrsm"))

    const Grid& g = RLPre.Grid();
    DistMatrix<F,STAR,STAR> RL(g);
    DistMatrix<F,VR,  STAR> RR(g);
    Copy( RLPre, RL, READ_PROXY );
    Copy( RRPre, RR, READ_WRITE_PROXY );

    PseudoTrsm( RL.Matrix(), RR.Matrix(), tol );

    Copy( RR, RRPre, RESTORE_READ_WRITE_PROXY );
}

// On output, the matrix Z contains the non-trivial portion of the interpolation
// matrix, and p contains the pivots used during the iterations of 
// pivoted QR. The input matrix A is unchanged.

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<Int>& p, 
  Matrix<F>& Z, const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("id::BusingerGolub"))
    typedef Base<F> Real;
    const Int n = A.Width();

    // Perform the pivoted QR factorization
    const Int numSteps = QR( A, p, ctrl );

    const Real eps = lapack::MachineEpsilon<Real>();
    const Real pinvTol = ( ctrl.adaptive ? ctrl.tol : numSteps*eps );

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    auto Rl = LockedView( A, IndexRange(0,numSteps), IndexRange(0,numSteps) );
    auto RR = LockedView( A, IndexRange(0,numStpes), IndexRange(numStpes,n) );
    Z = RR;
    PseudoTrsm( RL, Z, pinvTol );
}

template<typename F> 
inline void
BusingerGolub
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<Int>& p, 
  AbstractDistMatrix<F>& Z, const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("id::BusingerGolub"))
    typedef Base<F> Real;

    const Grid& g = APre.Grid();
    DistMatrix<F> A(g);
    Copy( APre, A, READ_WRITE_PROXY );

    const Int n = A.Width();

    // Perform the pivoted QR factorization
    const Int numSteps = QR( A, p, ctrl );

    const Real eps = lapack::MachineEpsilon<Real>();
    const Real pinvTol = ( ctrl.adaptive ? ctrl.tol : numSteps*eps );

    // Now form a minimizer of || RL Z - RR ||_2 via pseudo triangular solves
    auto RL = LockedViewRange( A, 0, 0,        numSteps, numSteps );
    auto RR = LockedViewRange( A, 0, numSteps, numSteps, n        );
    Copy( RR, Z );
    PseudoTrsm( RL, Z, pinvTol );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
}

} // namespace id

template<typename F> 
void ID
( const Matrix<F>& A, Matrix<Int>& p, 
        Matrix<F>& Z, const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    Matrix<F> B( A );
    id::BusingerGolub( B, p, Z, ctrl );
}

template<typename F> 
void ID
( Matrix<F>& A, Matrix<Int>& p, 
  Matrix<F>& Z, const QRCtrl<Base<F>> ctrl, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    id::BusingerGolub( B, p, Z, ctrl );
}

template<typename F> 
void ID
( const AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p, 
        AbstractDistMatrix<F>& Z, const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    DistMatrix<F> B( A );
    id::BusingerGolub( B, p, Z, ctrl );
}

template<typename F> 
void ID
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p, 
  AbstractDistMatrix<F>& Z, const QRCtrl<Base<F>> ctrl, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("ID"))
    if( canOverwrite )
    {
        id::BusingerGolub( A, p, Z, ctrl );
    }
    else
    {
        DistMatrix<F> B( A );
        id::BusingerGolub( B, p, Z, ctrl );
    }
}

#define PROTO(F) \
  template void ID \
  ( const Matrix<F>& A, Matrix<Int>& p, \
          Matrix<F>& Z, const QRCtrl<Base<F>> ctrl ); \
  template void ID \
  ( const AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p, \
          AbstractDistMatrix<F>& Z, const QRCtrl<Base<F>> ctrl ); \
  template void ID \
  ( Matrix<F>& A, Matrix<Int>& p, \
    Matrix<F>& Z, const QRCtrl<Base<F>> ctrl, bool canOverwrite ); \
  template void ID \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p, \
    AbstractDistMatrix<F>& Z, const QRCtrl<Base<F>> ctrl, bool canOverwrite ); 

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
