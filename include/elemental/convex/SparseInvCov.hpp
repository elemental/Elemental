/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SPARSEINVCOV_HPP
#define ELEM_SPARSEINVCOV_HPP

#include ELEM_AXPY_INC
#include ELEM_HILBERTSCHMIDT_INC
#include ELEM_SCALE_INC
#include ELEM_MAKEHERMITIAN_INC
#include ELEM_DETERMINANT_INC
#include ELEM_ENTRYWISEONENORM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_HERMITIANFUNCTION_INC
#include ELEM_COVARIANCE_INC
#include ELEM_SOFTTHRESHOLD_INC
#include ELEM_ZEROS_INC

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/covsel/covsel.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This ADMM attempts to solve the problem:
//     minimize Tr(S*X) - log det X + lambda ||X||_1
// where S is the empirical covariance of the data matrix D.
//

namespace elem {

template<typename F>
inline Int
SparseInvCov
( const Matrix<F>& D, Matrix<F>& X, Matrix<F>& Z, Matrix<F>& U,
  Base<F> lambda, Base<F> rho=1., Base<F> alpha=1.2, Int maxIter=500, 
  Base<F> absTol=1e-6, Base<F> relTol=1e-4, bool progress=true )
{
    DEBUG_ONLY(CallStackEntry cse("SparseInvCov"))
    typedef Base<F> Real;
    const Int n = D.Width();

    Matrix<F> S;
    Covariance( D, S );
    MakeHermitian( LOWER, S );

    Zeros( X, n, n );
    Zeros( Z, n, n );
    Zeros( U, n, n );
    
    Int numIter=0;
    Matrix<F> ZOld, XHat, T;
    while( numIter < maxIter )
    {
        ZOld = Z;

        // X := rho*(Z-U) - S
        X = Z;
        Axpy( F(-1), U, X );
        Scale( rho, X );
        Axpy( F(-1), S, X );

        // X := f(X), f(gamma) = (gamma+sqrt(gamma+4*rho)) / (2*rho)
        RealHermitianFunction
        ( LOWER, X, 
          [rho](Real gamma){return (gamma+Sqrt(gamma*gamma+4*rho))/(2*rho);} );
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
            const Real ZOne = EntrywiseOneNorm( Z );
            const Real objective = trace-safeDet.kappa*safeDet.n+lambda*ZOne;
            std::cout << numIter << ": "
              << "||X-Z||_F=" << rNorm << ", "
              << "epsPri=" << epsPri << ", "
              << "|rho| ||Z-ZOld||_F=" << sNorm << ", "
              << "epsDual=" << epsDual << ", "
              << "objective=" << objective << std::endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter )
        std::cout << "ADMM failed to converge" << std::endl;
    return numIter;
}

template<typename F>
inline Int
SparseInvCov
( const DistMatrix<F>& D, DistMatrix<F>& X, DistMatrix<F>& Z, DistMatrix<F>& U,
  Base<F> lambda, Base<F> rho=1., Base<F> alpha=1.2, Int maxIter=500, 
  Base<F> absTol=1e-6, Base<F> relTol=1e-4, bool progress=true )
{
    DEBUG_ONLY(CallStackEntry cse("SparseInvCov"))
    typedef Base<F> Real;
    const Grid& g = D.Grid();
    const Int n = D.Width();

    DistMatrix<F> S(g);
    Covariance( D, S );
    MakeHermitian( LOWER, S );

    Zeros( X, n, n );
    Zeros( Z, n, n );
    Zeros( U, n, n );
    
    Int numIter=0;
    DistMatrix<F> ZOld(g), XHat(g), T(g);
    while( numIter < maxIter )
    {
        ZOld = Z;

        // X := rho*(Z-U) - S
        X = Z;
        Axpy( F(-1), U, X );
        Scale( rho, X );
        Axpy( F(-1), S, X );

        // X := f(X), f(gamma) = (gamma+sqrt(gamma+4*rho)) / (2*rho)
        RealHermitianFunction
        ( LOWER, X, 
          [rho](Real gamma){return (gamma+Sqrt(gamma*gamma+4*rho))/(2*rho);} );
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
            const Real ZOne = EntrywiseOneNorm( Z );
            const Real objective = trace-safeDet.kappa*safeDet.n+lambda*ZOne;
            if( g.Rank() == 0 )
                std::cout << numIter << ": "
                  << "||X-Z||_F=" << rNorm << ", "
                  << "epsPri=" << epsPri << ", "
                  << "|rho| ||Z-ZOld||_F=" << sNorm << ", "
                  << "epsDual=" << epsDual << ", "
                  << "objective=" << objective << std::endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter && g.Rank() == 0 )
        std::cout << "ADMM failed to converge" << std::endl;
    return numIter;
}

} // namespace elem

#endif // ifndef ELEM_SPARSEINVCOV_HPP
