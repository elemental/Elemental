/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/basis_pursuit/basis_pursuit.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// Basis pursuit seeks the solution to A x = b which minimizes || x ||_1

namespace El {

template<typename F>
Int BasisPursuit
( const Matrix<F>& A, const Matrix<F>& b, Matrix<F>& z, 
  Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, 
  bool usePinv, Base<F> pinvTol, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("BasisPursuit"))
    // Find a means of quickly applyinv pinv(A) and then form pinv(A) b
    // NOTE: If m >= n and A has full column rank, then basis pursuit is 
    //       irrelevant, as there is a unique solution, which is found 
    //       through least squares. If A does *not* have full column rank,
    //       then the QR factorization is not enough.
    //       For the same reason, the LQ factorization will fail if A does
    //       not have full row rank.
    //
    // TODO: Instead form a basis for the null-space of A
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
    Matrix<F> x, u, t, zOld, xHat;
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
Int BasisPursuit
( const DistMatrix<F>& A, const DistMatrix<F>& b, DistMatrix<F>& z,
  Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, 
  bool usePinv, Base<F> pinvTol, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("BasisPursuit"))
    // Find a means of quickly applyinv pinv(A) and then form pinv(A) b
    // NOTE: If m >= n and A has full column rank, then basis pursuit is 
    //       irrelevant, as there is a unique solution, which is found 
    //       through least squares. If A does *not* have full column rank,
    //       then the QR factorization is not enough.
    //       For the same reason, the LQ factorization will fail if A does
    //       not have full row rank.
    //
    // TODO: Instead form a basis for the null-space of A
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
    DistMatrix<F> x(grid), u(grid), t(grid), zOld(grid), xHat(grid);
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

#define PROTO(F) \
  template Int BasisPursuit \
  ( const Matrix<F>& A, const Matrix<F>& b, Matrix<F>& z, \
    Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, \
    bool usePinv, Base<F> pinvTol, bool progress ); \
  template Int BasisPursuit \
  ( const DistMatrix<F>& A, const DistMatrix<F>& b, DistMatrix<F>& z, \
    Base<F> rho, Base<F> alpha, Int maxIter, Base<F> absTol, Base<F> relTol, \
    bool usePinv, Base<F> pinvTol, bool progress );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namepace elem
