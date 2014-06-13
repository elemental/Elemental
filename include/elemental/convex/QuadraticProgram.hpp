/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QUADRATICPROGRAM_HPP
#define ELEM_QUADRATICPROGRAM_HPP

#include ELEM_CLIP_INC
#include ELEM_DOT_INC
#include ELEM_MAKETRIANGULAR_INC
#include ELEM_UPDATEDIAGONAL_INC
#include ELEM_GEMV_INC
#include ELEM_HEMV_INC
#include ELEM_TRSV_INC
#include ELEM_CHOLESKY_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_TRIANGULARINVERSE_INC
#include ELEM_ZEROS_INC

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This ADMM attempts to solve the following quadratic program:
//     minimize    (1/2) x' P x + q' x
//     subject to  lb <= x <= ub
//

namespace elem {

template<typename Real>
inline Int
QuadraticProgram
( const Matrix<Real>& P, const Matrix<Real>& q, Real lb, Real ub,
  Matrix<Real>& x, Matrix<Real>& z, Matrix<Real>& u,
  Real rho=1., Real alpha=1.2, Int maxIter=500, 
  Real absTol=1e-6, Real relTol=1e-4, bool inv=false, bool progress=true )
{
    DEBUG_ONLY(CallStackEntry cse("QuadraticProgram"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    const Int n = P.Height();

    // Cache the factorization of P + rho*I
    Matrix<Real> LMod( P );
    UpdateDiagonal( LMod, rho );
    Cholesky( LOWER, LMod );
    MakeTriangular( LOWER, LMod );

    // Optionally invert the factor in place
    if( inv )
        TriangularInverse( LOWER, NON_UNIT, LMod );

    // Start the ADMM
    Int numIter=0;
    Matrix<Real> t, zOld, xHat;
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    Zeros( t, n, 1 );
    while( numIter < maxIter )
    {
        zOld = z;

        // x := (P+rho*I)^{-1} (rho(z-u)-q)
        x = z;
        Axpy( Real(-1), u, x );
        Scale( rho, x );
        Axpy( Real(-1), q, x );
        if( inv )
        {
            // TODO: Trmv
            Gemv( NORMAL, Real(1), LMod, x, t );
            Gemv( ADJOINT, Real(1), LMod, t, x );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, LMod, x );
            Trsv( LOWER, ADJOINT, NON_UNIT, LMod, x );
        }

        // xHat := alpha*x + (1-alpha)*zOld
        xHat = x;
        Scale( alpha, xHat );
        Axpy( 1-alpha, zOld, xHat );

        // z := Clip(xHat+u,lb,ub)
        z = xHat;
        Axpy( Real(1), u, z );
        Clip( z, lb, ub );

        // u := u + (xHat-z)
        Axpy( Real(1),  xHat, u );
        Axpy( Real(-1), z,    u );

        // Form (1/2) x' P x + q' x
        Zeros( t, n, 1 );
        Hemv( LOWER, Real(1), P, x, Real(0), t );
        const Real objective = Dot(x,t)/2 + Dot(q,x);

        // rNorm := || x - z ||_2
        t = x;
        Axpy( Real(-1), z, t );
        const Real rNorm = FrobeniusNorm( t );
        // sNorm := |rho| || z - zOld ||_2
        t = z;
        Axpy( Real(-1), zOld, t );
        const Real sNorm = Abs(rho)*FrobeniusNorm( t );

        const Real epsPri = Sqrt(Real(n))*absTol +
            relTol*Max(FrobeniusNorm(x),FrobeniusNorm(z));
        const Real epsDual = Sqrt(Real(n))*absTol +
            relTol*Abs(rho)*FrobeniusNorm(u);

        if( progress )
        {
            t = x;
            Clip( t, lb, ub );
            Axpy( Real(-1), x, t );
            const Real clipDist = FrobeniusNorm( t );
            std::cout << numIter << ": "
              << "||x-z||_2=" << rNorm << ", "
              << "epsPri=" << epsPri << ", "
              << "|rho| ||z-zOld||_2=" << sNorm << ", "
              << "epsDual=" << epsDual << ", "
              << "||x-Clip(x,lb,ub)||_2=" << clipDist << ", "
              << "(1/2) x' P x + q' x=" << objective << std::endl;
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
inline Int
QuadraticProgram
( const DistMatrix<Real>& P, const DistMatrix<Real>& q, Real lb, Real ub,
  DistMatrix<Real>& x, DistMatrix<Real>& z, DistMatrix<Real>& u,
  Real rho=1., Real alpha=1.2, Int maxIter=500, Real absTol=1e-6, 
  Real relTol=1e-4, bool inv=true, bool progress=true )
{
    DEBUG_ONLY(CallStackEntry cse("QuadraticProgram"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    const Grid& grid = P.Grid();
    const Int n = P.Height();

    // Cache the factorization of P + rho*I
    DistMatrix<Real> LMod( P );
    UpdateDiagonal( LMod, rho );
    Cholesky( LOWER, LMod );
    MakeTriangular( LOWER, LMod );

    // Optionally invert the factor in place
    if( inv )
        TriangularInverse( LOWER, NON_UNIT, LMod );

    // Start the ADMM
    Int numIter=0;
    DistMatrix<Real> t(grid), zOld(grid), xHat(grid);
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    Zeros( t, n, 1 );
    while( numIter < maxIter )
    {
        zOld = z;

        // x := (P+rho*I)^{-1} (rho(z-u)-q)
        x = z;
        Axpy( Real(-1), u, x );
        Scale( rho, x );
        Axpy( Real(-1), q, x );
        if( inv )
        {
            // TODO: Trmv
            Gemv( NORMAL, Real(1), LMod, x, t );
            Gemv( ADJOINT, Real(1), LMod, t, x );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, LMod, x );
            Trsv( LOWER, ADJOINT, NON_UNIT, LMod, x );
        }

        // xHat := alpha*x + (1-alpha)*zOld
        xHat = x;
        Scale( alpha, xHat );
        Axpy( 1-alpha, zOld, xHat );

        // z := Clip(xHat+u,lb,ub)
        z = xHat;
        Axpy( Real(1), u, z );
        Clip( z, lb, ub );

        // u := u + (xHat-z)
        Axpy( Real(1),  xHat, u );
        Axpy( Real(-1), z,    u );

        // Form (1/2) x' P x + q' x
        Zeros( t, n, 1 );
        Hemv( LOWER, Real(1), P, x, Real(0), t );
        const Real objective = Dot(x,t)/2 + Dot(q,x);

        // rNorm := || x - z ||_2
        t = x;
        Axpy( Real(-1), z, t );
        const Real rNorm = FrobeniusNorm( t );
        // sNorm := |rho| || z - zOld ||_2
        t = z;
        Axpy( Real(-1), zOld, t );
        const Real sNorm = Abs(rho)*FrobeniusNorm( t );

        const Real epsPri = Sqrt(Real(n))*absTol +
            relTol*Max(FrobeniusNorm(x),FrobeniusNorm(z));
        const Real epsDual = Sqrt(Real(n))*absTol +
            relTol*Abs(rho)*FrobeniusNorm(u);

        if( progress )
        {
            t = x;
            Clip( t, lb, ub );
            Axpy( Real(-1), x, t );
            const Real clipDist = FrobeniusNorm( t );
            if( grid.Rank() == 0 )
                std::cout << numIter << ": "
                  << "||x-z||_2=" << rNorm << ", "
                  << "epsPri=" << epsPri << ", "
                  << "|rho| ||z-zOld||_2=" << sNorm << ", "
                  << "epsDual=" << epsDual << ", "
                  << "||x-Clip(x,lb,ub)||_2=" << clipDist << ", "
                  << "(1/2) x' P x + q' x=" << objective << std::endl;
        }
        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter && grid.Rank() == 0 )
        std::cout << "ADMM failed to converge" << std::endl;
    return numIter;
}

} // namespace elem

#endif // ifndef ELEM_QUADRATICPROGRAM_HPP
