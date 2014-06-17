/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_IDENTITY_INC

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/lasso/lasso.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// The Least Absolute Shrinkage and Selection Operator (LASSO) solves the
// problem
//     min 1/2 || A x - b ||_2^2 + lambda || x ||_1
//
// Note that the prox map
//     min 1/2 || A x - b ||_2^2 + rho/2 || x - x_0 ||_2
//      x
// leads to the solution of the linear system
//     (A' A + rho) x = (A' b + rho x_0).
// When A is wider than it is tall, it is typically worthwhile to use the
// Woodbury matrix identity to re-express inv(A' A + rho) as 
//   (I - A' inv(A A' + rho) A) / rho.

namespace El {

template<typename F>
Int Lasso
( const Matrix<F>& A, const Matrix<F>& b, Base<F> lambda, Matrix<F>& z, 
  Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("Lasso"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<F> P;
    if( m >= n )
    {
        Identity( P, n, n );        
        Herk( LOWER, ADJOINT, F(1), A, F(rho), P );
    }
    else
    {
        Identity( P, m, m );
        Herk( LOWER, NORMAL, F(1), A, F(rho), P );
    }
    if( inv )
        HPDInverse( LOWER, P );
    else
        Cholesky( LOWER, P ); 

    // Cache w := A^H b
    Matrix<F> w;
    Gemv( ADJOINT, F(1), A, b, w );

    // Start the LASSO
    Int numIter=0;
    Matrix<F> x, u, s, zOld, xHat;
    Zeros( x, n, 1 );
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    while( numIter < maxIter )
    {
        zOld = z;

        // x := (A^H A + rho) \ (A^H b + rho*(z-u))
        x = w;
        Axpy(  rho, z, x );
        Axpy( -rho, u, x );
        if( m >= n )
        {
            if( inv )
            {
                s = x;
                Hemv( LOWER, F(1), P, s, F(0), x );
            }
            else
            {
                Trsv( LOWER, NORMAL, NON_UNIT, P, x );
                Trsv( LOWER, ADJOINT, NON_UNIT, P, x );
            }
        }
        else
        {
            Gemv( NORMAL, F(1), A, x, s );
            if( inv )
            {
                auto t( s );
                Hemv( LOWER, F(1), P, t, F(0), s );
            }
            else
            {
                Trsv( LOWER, NORMAL, NON_UNIT, P, s );
                Trsv( LOWER, ADJOINT, NON_UNIT, P, s );
            }
            Gemv( ADJOINT, F(-1), A, s, F(1), x );
            Scale( 1/rho, x );
        }

        // xHat := alpha x + (1-alpha) zOld
        xHat = x;
        Scale( alpha, xHat );
        Axpy( 1-alpha, zOld, xHat );

        // z := SoftThresh(xHat+u,lambda/rho)
        z = xHat;
        Axpy( F(1), u, z );
        SoftThreshold( z, lambda/rho );

        // u := u + (xHat - z)
        Axpy( F(1),  xHat, u );
        Axpy( F(-1), z,    u );

        // rNorm := || x - z ||_2
        s = x;
        Axpy( F(-1), z, s );
        const Real rNorm = FrobeniusNorm( s );

        // sNorm := || rho*(z-zOld) ||_2
        s = z;
        Axpy( F(-1), zOld, s );
        const Real sNorm = Abs(rho)*FrobeniusNorm( s );

        const Real epsPri = Sqrt(Real(n))*absTol +
            relTol*Max(FrobeniusNorm(x),FrobeniusNorm(z));
        const Real epsDual = Sqrt(Real(n))*absTol +
            relTol*Abs(rho)*FrobeniusNorm(u);

        if( progress )
        {
            s = b;
            Gemv( NORMAL, F(-1), A, x, F(1), s );
            const Real resid = FrobeniusNorm( s );
            const Real obj = Real(1)/Real(2)*resid*resid + lambda*OneNorm(z);
            std::cout << numIter << ": ||x-z||_2=" << rNorm
                      << ", epsPri=" << epsPri
                      << ", |rho| ||z-zOld||_2=" << sNorm
                      << ", and epsDual=" << epsDual << ", objective="
                      << obj << std::endl;
        }

        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter )
        std::cout << "Lasso failed to converge" << std::endl;
    return numIter;
}

template<typename F>
Int Lasso
( const DistMatrix<F>& A, const DistMatrix<F>& b, Base<F> lambda,
  DistMatrix<F>& z, 
  Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("Lasso"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();

    DistMatrix<F> P(g);
    if( m >= n )
    {
        Identity( P, n, n );        
        Herk( LOWER, ADJOINT, F(1), A, F(rho), P );
    }
    else
    {
        Identity( P, m, m );
        Herk( LOWER, NORMAL, F(1), A, F(rho), P );
    }
    if( inv )
        HPDInverse( LOWER, P );
    else
        Cholesky( LOWER, P ); 

    // Cache w := A^H b
    DistMatrix<F> w(g);
    Gemv( ADJOINT, F(1), A, b, w );

    // Start the LASSO
    Int numIter=0;
    DistMatrix<F> x(g), u(g), s(g), zOld(g), xHat(g);
    Zeros( x, n, 1 );
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    while( numIter < maxIter )
    {
        zOld = z;

        // x := (A^H A + rho) \ (A^H b + rho*(z-u))
        x = w;
        Axpy(  rho, z, x );
        Axpy( -rho, u, x );
        if( m >= n )
        {
            if( inv )
            {
                s = x;
                Hemv( LOWER, F(1), P, s, F(0), x );
            }
            else
            {
                Trsv( LOWER, NORMAL, NON_UNIT, P, x );
                Trsv( LOWER, ADJOINT, NON_UNIT, P, x );
            }
        }
        else
        {
            Gemv( NORMAL, F(1), A, x, s );
            if( inv )
            {
                auto t( s );
                Hemv( LOWER, F(1), P, t, F(0), s );
            }
            else
            {
                Trsv( LOWER, NORMAL, NON_UNIT, P, s );
                Trsv( LOWER, ADJOINT, NON_UNIT, P, s );
            }
            Gemv( ADJOINT, F(-1), A, s, F(1), x );
            Scale( 1/rho, x );
        }

        // xHat := alpha x + (1-alpha) zOld
        xHat = x;
        Scale( alpha, xHat );
        Axpy( 1-alpha, zOld, xHat );

        // z := SoftThresh(xHat+u,lambda/rho)
        z = xHat;
        Axpy( F(1), u, z );
        SoftThreshold( z, lambda/rho );

        // u := u + (xHat - z)
        Axpy( F(1),  xHat, u );
        Axpy( F(-1), z,    u );

        // rNorm := || x - z ||_2
        s = x;
        Axpy( F(-1), z, s );
        const Real rNorm = FrobeniusNorm( s );

        // sNorm := || rho*(z-zOld) ||_2
        s = z;
        Axpy( F(-1), zOld, s );
        const Real sNorm = Abs(rho)*FrobeniusNorm( s );

        const Real epsPri = Sqrt(Real(n))*absTol +
            relTol*Max(FrobeniusNorm(x),FrobeniusNorm(z));
        const Real epsDual = Sqrt(Real(n))*absTol +
            relTol*Abs(rho)*FrobeniusNorm(u);

        if( progress )
        {
            s = b;
            Gemv( NORMAL, F(-1), A, x, F(1), s );
            const Real resid = FrobeniusNorm( s );
            const Real obj = Real(1)/Real(2)*resid*resid + lambda*OneNorm(z);
            if( g.Rank() == 0 )
            {
                std::cout << numIter << ": ||x-z||_2=" << rNorm
                          << ", epsPri=" << epsPri
                          << ", |rho| ||z-zOld||_2=" << sNorm
                          << ", and epsDual=" << epsDual << ", objective="
                          << obj << std::endl;
            }
        }

        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter )
        std::cout << "Lasso failed to converge" << std::endl;
    return numIter;
}

#define PROTO(F) \
  template Int Lasso \
  ( const Matrix<F>& A, const Matrix<F>& b, Base<F> lambda, \
    Matrix<F>& z, \
    Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, \
    bool inv, bool progress ); \
  template Int Lasso \
  ( const DistMatrix<F>& A, const DistMatrix<F>& b, Base<F> lambda, \
    DistMatrix<F>& z, \
    Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, \
    bool inv, bool progress );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namepace elem
