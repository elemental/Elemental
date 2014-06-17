/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This ADMM attempts to solve the quadratic programs:
//     minimize    (1/2) x_j' P x_j + s_j' x_j
//     subject to  lb <= x_j <= ub
// where s_j is the j'th column of S and x_j is the corresponding solution 
// vector.
//

namespace El {

template<typename Real>
Int QuadraticProgram
( const Matrix<Real>& P, const Matrix<Real>& S, Real lb, Real ub,
  Matrix<Real>& Z, 
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("QuadraticProgram"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    const Int n = P.Height();
    const Int k = S.Width();

    // Cache the factorization of P + rho*I
    Matrix<Real> LMod( P );
    UpdateDiagonal( LMod, rho );
    if( inv )
    {
        HPDInverse( LOWER, LMod );
    }
    else
    {
        Cholesky( LOWER, LMod );
        MakeTriangular( LOWER, LMod );
    }

    // Start the ADMM
    Int numIter=0;
    Matrix<Real> X, U, T, ZOld, XHat;
    Zeros( Z, n, k );
    Zeros( U, n, k );
    Zeros( T, n, k );
    while( numIter < maxIter )
    {
        ZOld = Z;

        // x := (P+rho*I)^{-1} (rho(z-u)-q)
        X = Z;
        Axpy( Real(-1), U, X );
        Scale( rho, X );
        Axpy( Real(-1), S, X );
        if( inv )
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
        Scale( alpha, XHat );
        Axpy( 1-alpha, ZOld, XHat );

        // z := Clip(xHat+u,lb,ub)
        Z = XHat;
        Axpy( Real(1), U, Z );
        Clip( Z, lb, ub );

        // u := u + (xHat-z)
        Axpy( Real(1),  XHat, U );
        Axpy( Real(-1), Z,    U );

        // rNorm := || x - z ||_2
        T = X;
        Axpy( Real(-1), Z, T );
        const Real rNorm = FrobeniusNorm( T );

        // sNorm := |rho| || z - zOld ||_2
        T = Z;
        Axpy( Real(-1), ZOld, T );
        const Real sNorm = Abs(rho)*FrobeniusNorm( T );

        const Real epsPri = Sqrt(Real(n))*absTol +
            relTol*Max(FrobeniusNorm(X),FrobeniusNorm(Z));
        const Real epsDual = Sqrt(Real(n))*absTol +
            relTol*Abs(rho)*FrobeniusNorm(U);

        if( progress )
        {
            // Form (1/2) x' P x + s' x
            Zeros( T, n, k );
            Hemm( LEFT, LOWER, Real(1), P, X, Real(0), T );
            const Real objective = HilbertSchmidt(X,T)/2 + HilbertSchmidt(S,X);

            T = X;
            Clip( T, lb, ub );
            Axpy( Real(-1), X, T );
            const Real clipDist = FrobeniusNorm( T );
            std::cout << numIter << ": "
              << "||X-Z||_F=" << rNorm << ", "
              << "epsPri=" << epsPri << ", "
              << "|rho| ||Z-ZOld||_F=" << sNorm << ", "
              << "epsDual=" << epsDual << ", "
              << "||X-Clip(X,lb,ub)||_F=" << clipDist << ", "
              << "(1/2) <X,P X> + <C,X>=" << objective << std::endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter )
        std::cout << "ADMM failed to converge" << std::endl;
    return numIter;
}

template<typename Real>
Int QuadraticProgram
( const DistMatrix<Real>& P, const DistMatrix<Real>& S, Real lb, Real ub,
  DistMatrix<Real>& Z, 
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("QuadraticProgram"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    const Grid& grid = P.Grid();
    const Int n = P.Height();
    const Int k = S.Width();

    // Cache the factorization of P + rho*I
    DistMatrix<Real> LMod( P );
    UpdateDiagonal( LMod, rho );
    if( inv )
    {
        HPDInverse( LOWER, LMod );
    }
    else
    {
        Cholesky( LOWER, LMod );
        MakeTriangular( LOWER, LMod );
    }

    // Start the ADMM
    Int numIter=0;
    DistMatrix<Real> X(grid), U(grid), T(grid), ZOld(grid), XHat(grid);
    Zeros( Z, n, k );
    Zeros( U, n, k );
    Zeros( T, n, k );
    while( numIter < maxIter )
    {
        ZOld = Z;

        // x := (P+rho*I)^{-1} (rho(z-u)-q)
        X = Z;
        Axpy( Real(-1), U, X );
        Scale( rho, X );
        Axpy( Real(-1), S, X );
        if( inv )
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
        Scale( alpha, XHat );
        Axpy( 1-alpha, ZOld, XHat );

        // z := Clip(xHat+u,lb,ub)
        Z = XHat;
        Axpy( Real(1), U, Z );
        Clip( Z, lb, ub );

        // u := u + (xHat-z)
        Axpy( Real(1),  XHat, U );
        Axpy( Real(-1), Z,    U );

        // rNorm := || x - z ||_2
        T = X;
        Axpy( Real(-1), Z, T );
        const Real rNorm = FrobeniusNorm( T );

        // sNorm := |rho| || z - zOld ||_2
        T = Z;
        Axpy( Real(-1), ZOld, T );
        const Real sNorm = Abs(rho)*FrobeniusNorm( T );

        const Real epsPri = Sqrt(Real(n))*absTol +
            relTol*Max(FrobeniusNorm(X),FrobeniusNorm(Z));
        const Real epsDual = Sqrt(Real(n))*absTol +
            relTol*Abs(rho)*FrobeniusNorm(U);

        if( progress )
        {
            // Form (1/2) x' P x + s' x
            Zeros( T, n, k );
            Hemm( LEFT, LOWER, Real(1), P, X, Real(0), T );
            const Real objective = HilbertSchmidt(X,T)/2 + HilbertSchmidt(S,X);

            T = X;
            Clip( T, lb, ub );
            Axpy( Real(-1), X, T );
            const Real clipDist = FrobeniusNorm( T );
            if( grid.Rank() == 0 )
                std::cout << numIter << ": "
                  << "||X-Z||_F=" << rNorm << ", "
                  << "epsPri=" << epsPri << ", "
                  << "|rho| ||Z-ZOld||_F=" << sNorm << ", "
                  << "epsDual=" << epsDual << ", "
                  << "||X-Clip(X,lb,ub)||_2=" << clipDist << ", "
                  << "(1/2) <X,P X> + <S,X>=" << objective << std::endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter && grid.Rank() == 0 )
        std::cout << "ADMM failed to converge" << std::endl;
    return numIter;
}

#define PROTO(Real) \
  template Int QuadraticProgram \
  ( const Matrix<Real>& P, const Matrix<Real>& S, Real lb, Real ub, \
    Matrix<Real>& Z, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress ); \
  template Int QuadraticProgram \
  ( const DistMatrix<Real>& P, const DistMatrix<Real>& S, Real lb, Real ub, \
    DistMatrix<Real>& Z, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress );

PROTO(float)
PROTO(double)

} // namespace El
