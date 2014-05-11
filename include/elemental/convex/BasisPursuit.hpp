/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BASISPURSUIT_HPP
#define ELEM_BASISPURSUIT_HPP

#include ELEM_GEMV_INC
#include ELEM_TRSV_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_LQ_INC
#include ELEM_QR_INC
#include ELEM_PSEUDOINVERSE_INC
#include ELEM_SOFTTHRESHOLD_INC
#include ELEM_ZEROS_INC

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/basis_pursuit/basis_pursuit.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// Basis pursuit seeks the solution to A x = b which minimizes || x ||_1

namespace elem {

template<typename F>
inline Int
BasisPursuit
( const Matrix<F>& A, const Matrix<F>& b,
  Matrix<F>& x, Matrix<F>& z, Matrix<F>& u, Base<F> rho=1., Base<F> alpha=1.2, 
  Int maxIter=500, Base<F> absTol=1e-6, Base<F> relTol=1e-4, bool usePinv=false,
  Base<F> pinvTol=0, bool progress=true )
{
    DEBUG_ONLY(CallStackEntry cse("BasisPursuit"))
    // Find a means of quickly applyinv pinv(A) and then form pinv(A) b
    // NOTE: If m >= n and A has full column rank, then basis pursuit is 
    //       irrelevant, as there is a unique solution, which is found 
    //       through least squares. If A does *not* have full column rank,
    //       then the QR factorization is not enough.
    //       For the same reason, the LQ factorization will fail if A does
    //       not have full row rank.
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<F> pinvA;
    Matrix<F> Q, L, R;
    Matrix<F> q, s;
    if( usePinv )
    {
        pinvA = A;
        Pseudoinverse( pinvA, pinvTol );
        Gemv( NORMAL, F(1), pinvA, b, q );
    }
    else if( m >= n )
    {
        Q = A;
        qr::Explicit( Q, R );
        Gemv( ADJOINT, F(1), Q, b, q );
        Trsv( UPPER, NORMAL, NON_UNIT, R, q );
    }
    else
    {
        Q = A;
        lq::Explicit( L, Q );
        s = b;
        Trsv( LOWER, NORMAL, NON_UNIT, L, s );
        Gemv( ADJOINT, F(1), Q, s, q );
    }

    if( progress )
    {
        const Real qOneNorm = OneNorm( q );
        std::cout << " || pinv(A) b ||_1 = " << qOneNorm << std::endl;
    }

    // Start the basis pursuit
    Int numIter=0;
    Matrix<F> t, zOld, xHat;
    Zeros( x, n, 1 );
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    while( numIter < maxIter )
    {
        zOld = z;

        // x := P*(z-u) + q
        //    = (I-pinv(A)*A)(z-u) + q
        //    = (z-u) - pinv(A)*A*(z-u) + q
        s = z;
        Axpy( F(-1), u, s );
        x = s;
        Gemv( NORMAL, F(1), A, s, t );
        if( usePinv )
        {
            Gemv( NORMAL, F(1), pinvA, t, s );
        }
        else if( m >= n )
        {
            Gemv( ADJOINT, F(1), Q, t, s );
            Trsv( UPPER, NORMAL, NON_UNIT, R, s );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, L, t );
            Gemv( ADJOINT, F(1), Q, t, s );
        }
        Axpy( F(-1), s, x );
        Axpy( F(1),  q, x );

        // xHat := alpha x + (1-alpha) zOld
        xHat = x;
        Scale( alpha, xHat );
        Axpy( 1-alpha, zOld, xHat );

        // z := SoftThresh(xHat+u,1/rho)
        z = xHat;
        Axpy( F(1), u, z );
        SoftThreshold( z, 1/rho );

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
            const Real xOneNorm = OneNorm( x );
            std::cout << numIter << ": ||x-z||_2=" << rNorm
                      << ", epsPri=" << epsPri
                      << ", |rho| ||z-zOld||_2=" << sNorm
                      << ", and epsDual=" << epsDual << ", ||x||_1="
                      << xOneNorm << std::endl;
        }

        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter )
        std::cout << "Basis pursuit failed to converge" << std::endl;
    return numIter;
}

template<typename F>
inline Int
BasisPursuit
( const DistMatrix<F>& A, const DistMatrix<F>& b, 
  DistMatrix<F>& x, DistMatrix<F>& z, DistMatrix<F>& u, 
  Base<F> rho=1., Base<F> alpha=1.2, 
  Int maxIter=500, Base<F> absTol=1e-6, Base<F> relTol=1e-4, bool usePinv=false,
  Base<F> pinvTol=0, bool progress=true )
{
    DEBUG_ONLY(CallStackEntry cse("BasisPursuit"))
    // Find a means of quickly applyinv pinv(A) and then form pinv(A) b
    // NOTE: If m >= n and A has full column rank, then basis pursuit is 
    //       irrelevant, as there is a unique solution, which is found 
    //       through least squares. If A does *not* have full column rank,
    //       then the QR factorization is not enough.
    //       For the same reason, the LQ factorization will fail if A does
    //       not have full row rank.
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    DistMatrix<F> pinvA(grid);
    DistMatrix<F> Q(grid), L(grid), R(grid);
    DistMatrix<F> q(grid), s(grid);
    if( usePinv )
    {
        pinvA = A;
        Pseudoinverse( pinvA, pinvTol );
        Gemv( NORMAL, F(1), pinvA, b, q );
    }
    else if( m >= n )
    {
        Q = A;
        qr::Explicit( Q, R );
        Gemv( ADJOINT, F(1), Q, b, q );
        Trsv( UPPER, NORMAL, NON_UNIT, R, q );
    }
    else
    {
        Q = A;
        lq::Explicit( L, Q );
        s = b;
        Trsv( LOWER, NORMAL, NON_UNIT, L, s );
        Gemv( ADJOINT, F(1), Q, s, q );
    }

    if( progress )
    {
        const Real qOneNorm = OneNorm( q );
        if( grid.Rank() == 0 )
            std::cout << " || pinv(A) b ||_1 = " << qOneNorm << std::endl;
    }

    // Start the basis pursuit
    Int numIter=0;
    DistMatrix<F> t(grid), zOld(grid), xHat(grid);
    Zeros( x, n, 1 );
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    while( numIter < maxIter )
    {
        zOld = z;

        // x := P*(z-u) + q
        //    = (I-pinv(A)*A)(z-u) + q
        //    = (z-u) - pinv(A)*A*(z-u) + q
        s = z;
        Axpy( F(-1), u, s );
        x = s;
        Gemv( NORMAL, F(1), A, s, t );
        if( usePinv )
        {
            Gemv( NORMAL, F(1), pinvA, t, s );
        }
        else if( m >= n )
        {
            Gemv( ADJOINT, F(1), Q, t, s );
            Trsv( UPPER, NORMAL, NON_UNIT, R, s );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, L, t );
            Gemv( ADJOINT, F(1), Q, t, s );
        }
        Axpy( F(-1), s, x );
        Axpy( F(1),  q, x );

        // xHat := alpha x + (1-alpha) zOld
        xHat = x;
        Scale( alpha, xHat );
        Axpy( 1-alpha, zOld, xHat );

        // z := SoftThresh(xHat+u,1/rho)
        z = xHat;
        Axpy( F(1), u, z );
        SoftThreshold( z, 1/rho );

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
            const Real xOneNorm = OneNorm( x );
            if( grid.Rank() == 0 )
                std::cout << numIter << ": ||x-z||_2=" << rNorm
                          << ", epsPri=" << epsPri
                          << ", |rho| ||z-zOld||_2=" << sNorm
                          << ", and epsDual=" << epsDual << ", ||x||_1="
                          << xOneNorm << std::endl;
        }

        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( maxIter == numIter && grid.Rank() == 0 )
        std::cout << "Basis pursuit failed to converge" << std::endl;
    return numIter;
}

} // namepace elem

#endif // ifndef ELEM_BASISPURSUIT_HPP
