/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// TODO: Add a conic-form ADMM (i.e., x >= 0)

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This ADMM attempts to solve the quadratic programs:
//     minimize    (1/2) x' Q x + c' x
//     subject to  lb <= x <= ub
// where c and x are corresponding columns of the matrices C and X, 
// respectively.
//

namespace El {
namespace qp {
namespace box {

template<typename Real,typename>
Int
ADMM
( const Matrix<Real>& Q,
  const Matrix<Real>& C, 
        Real lb,
        Real ub,
        Matrix<Real>& Z, 
  const ADMMCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int n = Q.Height();
    const Int k = C.Width();

    // Cache the factorization of Q + rho*I
    Matrix<Real> LMod( Q );
    ShiftDiagonal( LMod, ctrl.rho );
    if( ctrl.inv )
    {
        HPDInverse( LOWER, LMod );
    }
    else
    {
        Cholesky( LOWER, LMod );
        MakeTrapezoidal( LOWER, LMod );
    }

    // Start the ADMM
    Int numIter=0;
    Matrix<Real> X, U, T, ZOld, XHat;
    Zeros( Z, n, k );
    Zeros( U, n, k );
    Zeros( T, n, k );
    while( numIter < ctrl.maxIter )
    {
        ZOld = Z;

        // x := (Q+rho*I)^{-1} (rho(z-u)-q)
        X = Z;
        X -= U;
        X *= ctrl.rho;
        X -= C;
        if( ctrl.inv )
        {
            auto Y( X );
            Hemm( LEFT, LOWER, Real(1), LMod, Y, Real(0), X );
        }
        else
        {
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, Real(1), LMod, X );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, Real(1), LMod, X );
        }

        // xHat := alpha*x + (1-alpha)*zOld
        XHat = X;
        XHat *= ctrl.alpha;
        Axpy( 1-ctrl.alpha, ZOld, XHat );

        // z := Clip(xHat+u,lb,ub)
        Z = XHat;
        Z += U;
        Clip( Z, lb, ub );

        // u := u + (xHat-z)
        U += XHat;
        U -= Z;

        // rNorm := || x - z ||_2
        T = X;
        T -= Z;
        const Real rNorm = FrobeniusNorm( T );

        // sNorm := |rho| || z - zOld ||_2
        T = Z;
        T -= ZOld;
        const Real sNorm = Abs(ctrl.rho)*FrobeniusNorm( T );

        const Real epsPri = Sqrt(Real(n))*ctrl.absTol +
            ctrl.relTol*Max(FrobeniusNorm(X),FrobeniusNorm(Z));
        const Real epsDual = Sqrt(Real(n))*ctrl.absTol +
            ctrl.relTol*Abs(ctrl.rho)*FrobeniusNorm(U);

        if( ctrl.print )
        {
            // Form (1/2) x' Q x + c' x
            Zeros( T, n, k );
            Hemm( LEFT, LOWER, Real(1), Q, X, Real(0), T );
            const Real objective = HilbertSchmidt(X,T)/2 + HilbertSchmidt(C,X);

            T = X;
            Clip( T, lb, ub );
            T -= X;
            const Real clipDist = FrobeniusNorm( T );
            cout << numIter << ": "
              << "||X-Z||_F=" << rNorm << ", "
              << "epsPri=" << epsPri << ", "
              << "|rho| ||Z-ZOld||_F=" << sNorm << ", "
              << "epsDual=" << epsDual << ", "
              << "||X-Clip(X,lb,ub)||_F=" << clipDist << ", "
              << "(1/2) <X,Q X> + <C,X>=" << objective << endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( ctrl.maxIter == numIter )
        Output("ADMM failed to converge");
    return numIter;
}

template<typename Real,typename>
Int
ADMM
( const ElementalMatrix<Real>& QPre,
  const ElementalMatrix<Real>& CPre, 
        Real lb,
        Real ub,
        ElementalMatrix<Real>& ZPre, 
  const ADMMCtrl<Real>& ctrl )
{
    DEBUG_CSE

    DistMatrixReadProxy<Real,Real,MC,MR>
      QProx( QPre ),
      CProx( CPre );
    DistMatrixWriteProxy<Real,Real,MC,MR>
      ZProx( ZPre );
    auto& Q = QProx.GetLocked();
    auto& C = CProx.GetLocked();
    auto& Z = ZProx.Get();

    const Grid& grid = Q.Grid();
    const Int n = Q.Height();
    const Int k = C.Width();

    // Cache the factorization of Q + rho*I
    DistMatrix<Real> LMod( Q );
    ShiftDiagonal( LMod, ctrl.rho );
    if( ctrl.inv )
    {
        HPDInverse( LOWER, LMod );
    }
    else
    {
        Cholesky( LOWER, LMod );
        MakeTrapezoidal( LOWER, LMod );
    }

    // Start the ADMM
    Int numIter=0;
    DistMatrix<Real> X(grid), U(grid), T(grid), ZOld(grid), XHat(grid);
    Zeros( Z, n, k );
    Zeros( U, n, k );
    Zeros( T, n, k );
    while( numIter < ctrl.maxIter )
    {
        ZOld = Z;

        // x := (Q+rho*I)^{-1} (rho(z-u)-q)
        X = Z;
        X -= U;
        X *= ctrl.rho;
        X -= C;
        if( ctrl.inv )
        {
            auto Y( X );
            Hemm( LEFT, LOWER, Real(1), LMod, Y, Real(0), X );
        }
        else
        {
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, Real(1), LMod, X );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, Real(1), LMod, X );
        }

        // xHat := alpha*x + (1-alpha)*zOld
        XHat = X;
        XHat *= ctrl.alpha;
        Axpy( 1-ctrl.alpha, ZOld, XHat );

        // z := Clip(xHat+u,lb,ub)
        Z = XHat;
        Z += U;
        Clip( Z, lb, ub );

        // u := u + (xHat-z)
        U += XHat;
        U -= Z;

        // rNorm := || x - z ||_2
        T = X;
        T -= Z;
        const Real rNorm = FrobeniusNorm( T );

        // sNorm := |rho| || z - zOld ||_2
        T = Z;
        T -= ZOld;
        const Real sNorm = Abs(ctrl.rho)*FrobeniusNorm( T );

        const Real epsPri = Sqrt(Real(n))*ctrl.absTol +
            ctrl.relTol*Max(FrobeniusNorm(X),FrobeniusNorm(Z));
        const Real epsDual = Sqrt(Real(n))*ctrl.absTol +
            ctrl.relTol*Abs(ctrl.rho)*FrobeniusNorm(U);

        if( ctrl.print )
        {
            // Form (1/2) x' Q x + c' x
            Zeros( T, n, k );
            Hemm( LEFT, LOWER, Real(1), Q, X, Real(0), T );
            const Real objective = HilbertSchmidt(X,T)/2 + HilbertSchmidt(C,X);

            T = X;
            Clip( T, lb, ub );
            T -= X;
            const Real clipDist = FrobeniusNorm( T );
            if( grid.Rank() == 0 )
                cout << numIter << ": "
                  << "||X-Z||_F=" << rNorm << ", "
                  << "epsPri=" << epsPri << ", "
                  << "|rho| ||Z-ZOld||_F=" << sNorm << ", "
                  << "epsDual=" << epsDual << ", "
                  << "||X-Clip(X,lb,ub)||_2=" << clipDist << ", "
                  << "(1/2) <X,Q X> + <C,X>=" << objective << endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( ctrl.maxIter == numIter && grid.Rank() == 0 )
        Output("ADMM failed to converge");
    return numIter;
}

#define PROTO(Real) \
  template Int ADMM \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& C, \
          Real lb, \
          Real ub, \
          Matrix<Real>& Z, \
    const ADMMCtrl<Real>& ctrl ); \
  template Int ADMM \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& C, \
          Real lb, \
          Real ub, \
          ElementalMatrix<Real>& Z, \
    const ADMMCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace box
} // namespace qp
} // namespace El
