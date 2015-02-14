/*
   Copyright (c) 2009-2015, Jack Poulson
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
( const Matrix<F>& D, Base<F> lambda, Matrix<F>& Z,
  Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, 
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SparseInvCov"))
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
    while( numIter < maxIter )
    {
        ZOld = Z;

        // X := rho*(Z-U) - S
        X = Z;
        Axpy( F(-1), U, X );
        Scale( rho, X );
        Axpy( F(-1), S, X );

        // X := f(X), f(gamma) = (gamma+sqrt(gamma+4*rho)) / (2*rho)
        auto eigMap = 
          [rho](Real gamma){return (gamma+Sqrt(gamma*gamma+4*rho))/(2*rho);};
        HermitianFunction( LOWER, X, function<Real(Real)>(eigMap) );
        // Make X explicitly Hermitian since HermitianHilbertSchmidt is not
        // yet available. This should result in Z and U remaining explicitly
        // Hermitian.
        MakeHermitian( LOWER, X );

        // XHat := alpha*X + (1-alpha)*ZOld
        XHat = X;
        Scale( alpha, XHat );
        Axpy( 1-alpha, ZOld, XHat );

        // Z := SoftThreshold(XHat+U,lambda/rho)
        Z = XHat;
        Axpy( Real(1), U, Z );
        SoftThreshold( Z, lambda/rho );

        // U := U + (XHat-Z)
        Axpy( Real(1),  XHat, U );
        Axpy( Real(-1), Z,    U );

        // rNorm := || X - Z ||_F
        T = X;
        Axpy( Real(-1), Z, T );
        const Real rNorm = FrobeniusNorm(T);
        // sNorm := |rho| || Z - ZOld ||_F
        T = Z;
        Axpy( Real(-1), ZOld, T );
        const Real sNorm = Abs(rho)*FrobeniusNorm(T);

        const Real epsPri = n*absTol + 
            relTol*Max(FrobeniusNorm(X),FrobeniusNorm(Z));
        const Real epsDual = n*absTol + relTol*Abs(rho)*FrobeniusNorm(U);

        if( progress )
        {
            const Real trace = RealPart(HilbertSchmidt( S, X ));
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
    if( maxIter == numIter )
        cout << "ADMM failed to converge" << endl;
    return numIter;
}

template<typename F>
Int SparseInvCov
( const AbstractDistMatrix<F>& D, Base<F> lambda, AbstractDistMatrix<F>& ZPre,
  Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, 
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SparseInvCov"))

    auto ZPtr = WriteProxy<F,MC,MR>( &ZPre );
    auto& Z = *ZPtr;

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
    while( numIter < maxIter )
    {
        ZOld = Z;

        // X := rho*(Z-U) - S
        X = Z;
        Axpy( F(-1), U, X );
        Scale( rho, X );
        Axpy( F(-1), S, X );

        // X := f(X), f(gamma) = (gamma+sqrt(gamma+4*rho)) / (2*rho)
        auto eigMap = 
          [rho](Real gamma){return (gamma+Sqrt(gamma*gamma+4*rho))/(2*rho);};
        HermitianFunction( LOWER, X, function<Real(Real)>(eigMap) );
        // Make X explicitly Hermitian since HermitianHilbertSchmidt is not
        // yet available. This should result in Z and U remaining explicitly
        // Hermitian.
        MakeHermitian( LOWER, X );

        // XHat := alpha*X + (1-alpha)*ZOld
        XHat = X;
        Scale( alpha, XHat );
        Axpy( 1-alpha, ZOld, XHat );

        // Z := SoftThreshold(XHat+U,lambda/rho)
        Z = XHat;
        Axpy( Real(1), U, Z );
        SoftThreshold( Z, lambda/rho );

        // U := U + (XHat-Z)
        Axpy( Real(1),  XHat, U );
        Axpy( Real(-1), Z,    U );

        // rNorm := || X - Z ||_F
        T = X;
        Axpy( Real(-1), Z, T );
        const Real rNorm = FrobeniusNorm(T);
        // sNorm := |rho| || Z - ZOld ||_F
        T = Z;
        Axpy( Real(-1), ZOld, T );
        const Real sNorm = Abs(rho)*FrobeniusNorm(T);

        const Real epsPri = n*absTol + 
            relTol*Max(FrobeniusNorm(X),FrobeniusNorm(Z));
        const Real epsDual = n*absTol + relTol*Abs(rho)*FrobeniusNorm(U);

        if( progress )
        {
            const Real trace = RealPart(HilbertSchmidt( S, X ));
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
    if( maxIter == numIter && g.Rank() == 0 )
        cout << "ADMM failed to converge" << endl;
    return numIter;
}

#define PROTO(F) \
  template Int SparseInvCov \
  ( const Matrix<F>& D, Base<F> lambda, Matrix<F>& Z, \
    Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, \
    bool progress ); \
  template Int SparseInvCov \
  ( const AbstractDistMatrix<F>& D, Base<F> lambda, AbstractDistMatrix<F>& Z, \
    Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, \
    bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
