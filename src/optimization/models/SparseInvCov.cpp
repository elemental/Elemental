/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/covsel/covsel.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This ADMM attempts to solve the problem:
//     minimize Tr(S*X) - log det X + lambda ||X||_1
// where S is the empirical covariance of the data matrix D.
//

namespace El {

template<typename Field>
Int SparseInvCov
( const Matrix<Field>& D,
        Base<Field> lambda,
        Matrix<Field>& Z,
  const SparseInvCovCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = D.Width();

    Matrix<Field> S;
    Covariance( D, S );
    MakeHermitian( LOWER, S );

    Int numIter=0;
    Matrix<Field> X, U, ZOld, XHat, T;
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
          [&](const Real& gamma)
          { return (gamma+Sqrt(gamma*gamma+4*ctrl.rho))/(2*ctrl.rho); };
        HermitianFunction( LOWER, X, MakeFunction(eigMap) );
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
            Output
            (numIter,": ","||X-Z||_F=",rNorm,", epsPri=",epsPri,
             ", |rho| ||Z-ZOld||_F=",sNorm,", epsDual=",epsDual,
             ", objective=",objective);
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( ctrl.maxIter == numIter )
        RuntimeError("ADMM failed to converge");
    return numIter;
}

template<typename Field>
Int SparseInvCov
( const AbstractDistMatrix<Field>& D,
        Base<Field> lambda,
        AbstractDistMatrix<Field>& ZPre,
  const SparseInvCovCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixWriteProxy<Field,Field,MC,MR> ZProx( ZPre );
    auto& Z = ZProx.Get();

    typedef Base<Field> Real;
    const Grid& g = D.Grid();
    const Int n = D.Width();

    DistMatrix<Field> S(g);
    Covariance( D, S );
    MakeHermitian( LOWER, S );

    Int numIter=0;
    DistMatrix<Field> X(g), U(g), ZOld(g), XHat(g), T(g);
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
          [&](const Real& gamma)
          { return (gamma+Sqrt(gamma*gamma+4*ctrl.rho))/(2*ctrl.rho); };
        HermitianFunction( LOWER, X, MakeFunction(eigMap) );
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
                Output
                (numIter,": ","||X-Z||_F=",rNorm,", epsPri=",epsPri,
                 ", |rho| ||Z-ZOld||_F=",sNorm,", epsDual=",epsDual,
                 ", objective=",objective);
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( ctrl.maxIter == numIter )
        RuntimeError("ADMM failed to converge");
    return numIter;
}

#define PROTO(Field) \
  template Int SparseInvCov \
  ( const Matrix<Field>& D, \
          Base<Field> lambda, \
          Matrix<Field>& Z, \
    const SparseInvCovCtrl<Base<Field>>& ctrl ); \
  template Int SparseInvCov \
  ( const AbstractDistMatrix<Field>& D, \
          Base<Field> lambda, \
          AbstractDistMatrix<Field>& Z, \
    const SparseInvCovCtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#include <El/macros/Instantiate.h>

} // namespace El
