/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/covsel/covsel.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This ADMM attempts to solve the problem:
//     minimize Tr(S*X) - log det X + lambda ||X||_1
// where S is the empirical covariance of the data matrix D.
//

namespace El {

template<typename F>
Int SparseInvCov
( const Matrix<F>& D,
        Base<F> lambda,
        Matrix<F>& Z,
  const SparseInvCovCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SparseInvCov"))
    typedef Base<F> Real;
    const Int n = D.Width();

    Matrix<F> S;
    Covariance( D, S );
    MakeHermitian( LOWER, S );
   
    Int numIter=0;
    Matrix<F> X, U, ZOld, XHat, T;
    Zeros( X, n, n );
    Zeros( Z, n, n );
    Zeros( U, n, n );
    while( numIter < ctrl.maxIter )
    {
        ZOld = Z;

        // X := rho*(Z-U) - S
        X = Z;
        X -= U;
        X *= ctrl.rho;
        X -= S;

        // X := f(X), f(gamma) = (gamma+sqrt(gamma+4*rho)) / (2*rho)
        auto eigMap = 
          [ctrl](Real gamma)
          { return (gamma+Sqrt(gamma*gamma+4*ctrl.rho))/(2*ctrl.rho); };
        HermitianFunction( LOWER, X, function<Real(Real)>(eigMap) );
        // Make X explicitly Hermitian since HermitianHilbertSchmidt is not
        // yet available. This should result in Z and U remaining explicitly
        // Hermitian.
        MakeHermitian( LOWER, X );

        // XHat := alpha*X + (1-alpha)*ZOld
        XHat = X;
        XHat *= ctrl.alpha;
        Axpy( 1-ctrl.alpha, ZOld, XHat );

        // Z := SoftThreshold(XHat+U,lambda/rho)
        Z = XHat;
        Z += U;
        SoftThreshold( Z, lambda/ctrl.rho );

        // U := U + (XHat-Z)
        U += XHat;
        U -= Z;

        // rNorm := || X - Z ||_F
        T = X;
        T -= Z;
        const Real rNorm = FrobeniusNorm(T);
        // sNorm := |rho| || Z - ZOld ||_F
        T = Z;
        T -= ZOld;
        const Real sNorm = Abs(ctrl.rho)*FrobeniusNorm(T);

        const Real epsPri = n*ctrl.absTol + 
            ctrl.relTol*Max(FrobeniusNorm(X),FrobeniusNorm(Z));
        const Real epsDual = n*ctrl.absTol + 
            ctrl.relTol*Abs(ctrl.rho)*FrobeniusNorm(U);

        if( ctrl.progress )
        {
            const Real trace = RealPart(HilbertSchmidt(S,X));
            const SafeProduct<Real> safeDet = SafeHPDDeterminant( LOWER, X );
            const Real ZOne = EntrywiseNorm( Z, Real(1) );
            const Real objective = trace-safeDet.kappa*safeDet.n+lambda*ZOne;
            cout << numIter << ": "
              << "||X-Z||_F=" << rNorm << ", "
              << "epsPri=" << epsPri << ", "
              << "|rho| ||Z-ZOld||_F=" << sNorm << ", "
              << "epsDual=" << epsDual << ", "
              << "objective=" << objective << endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( ctrl.maxIter == numIter )
        cout << "ADMM failed to converge" << endl;
    return numIter;
}

template<typename F>
Int SparseInvCov
( const ElementalMatrix<F>& D,
        Base<F> lambda,
        ElementalMatrix<F>& ZPre,
  const SparseInvCovCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SparseInvCov"))

    DistMatrixWriteProxy<F,F,MC,MR> ZProx( ZPre );
    auto& Z = ZProx.Get();

    typedef Base<F> Real;
    const Grid& g = D.Grid();
    const Int n = D.Width();

    DistMatrix<F> S(g);
    Covariance( D, S );
    MakeHermitian( LOWER, S );
   
    Int numIter=0;
    DistMatrix<F> X(g), U(g), ZOld(g), XHat(g), T(g);
    Zeros( X, n, n );
    Zeros( Z, n, n );
    Zeros( U, n, n );
    while( numIter < ctrl.maxIter )
    {
        ZOld = Z;

        // X := rho*(Z-U) - S
        X = Z;
        X -= U;
        X *= ctrl.rho;
        X -= S;

        // X := f(X), f(gamma) = (gamma+sqrt(gamma+4*rho)) / (2*rho)
        auto eigMap = 
          [ctrl](Real gamma)
          { return (gamma+Sqrt(gamma*gamma+4*ctrl.rho))/(2*ctrl.rho); };
        HermitianFunction( LOWER, X, function<Real(Real)>(eigMap) );
        // Make X explicitly Hermitian since HermitianHilbertSchmidt is not
        // yet available. This should result in Z and U remaining explicitly
        // Hermitian.
        MakeHermitian( LOWER, X );

        // XHat := alpha*X + (1-alpha)*ZOld
        XHat = X;
        XHat *= ctrl.alpha;
        Axpy( 1-ctrl.alpha, ZOld, XHat );

        // Z := SoftThreshold(XHat+U,lambda/rho)
        Z = XHat;
        Z += U;
        SoftThreshold( Z, lambda/ctrl.rho );

        // U := U + (XHat-Z)
        U += XHat;
        U -= Z;

        // rNorm := || X - Z ||_F
        T = X;
        T -= Z;
        const Real rNorm = FrobeniusNorm(T);
        // sNorm := |rho| || Z - ZOld ||_F
        T = Z;
        T -= ZOld;
        const Real sNorm = Abs(ctrl.rho)*FrobeniusNorm(T);

        const Real epsPri = n*ctrl.absTol + 
            ctrl.relTol*Max(FrobeniusNorm(X),FrobeniusNorm(Z));
        const Real epsDual = n*ctrl.absTol + 
            ctrl.relTol*Abs(ctrl.rho)*FrobeniusNorm(U);

        if( ctrl.progress )
        {
            const Real trace = RealPart(HilbertSchmidt(S,X));
            const SafeProduct<Real> safeDet = SafeHPDDeterminant( LOWER, X );
            const Real ZOne = EntrywiseNorm( Z, Real(1) );
            const Real objective = trace-safeDet.kappa*safeDet.n+lambda*ZOne;
            if( g.Rank() == 0 )
                cout << numIter << ": "
                  << "||X-Z||_F=" << rNorm << ", "
                  << "epsPri=" << epsPri << ", "
                  << "|rho| ||Z-ZOld||_F=" << sNorm << ", "
                  << "epsDual=" << epsDual << ", "
                  << "objective=" << objective << endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( ctrl.maxIter == numIter && g.Rank() == 0 )
        cout << "ADMM failed to converge" << endl;
    return numIter;
}

#define PROTO(F) \
  template Int SparseInvCov \
  ( const Matrix<F>& D, \
          Base<F> lambda, \
          Matrix<F>& Z, \
    const SparseInvCovCtrl<Base<F>>& ctrl ); \
  template Int SparseInvCov \
  ( const ElementalMatrix<F>& D, \
          Base<F> lambda, \
          ElementalMatrix<F>& Z, \
    const SparseInvCovCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
