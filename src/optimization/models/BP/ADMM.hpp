/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// These implementations are adaptations of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/basis_pursuit/basis_pursuit.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// Basis pursuit seeks the solution to A x = b which minimizes || x ||_1

namespace El {
namespace bp {

template<typename Field>
Int ADMM
( const Matrix<Field>& A,
  const Matrix<Field>& b,
        Matrix<Field>& z,
  const ADMMCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    // Find a means of quickly applyinv pinv(A) and then form pinv(A) b
    // NOTE: If m >= n and A has full column rank, then basis pursuit is
    //       irrelevant, as there is a unique solution, which is found
    //       through least squares. If A does *not* have full column rank,
    //       then the QR factorization is not enough.
    //       For the same reason, the LQ factorization will fail if A does
    //       not have full row rank.
    //
    // TODO(poulson): Instead form a basis for the null-space of A
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Field> pinvA;
    Matrix<Field> Q, L, R;
    Matrix<Field> q, s;
    if( ctrl.usePinv )
    {
        pinvA = A;
        Pseudoinverse( pinvA, ctrl.pinvTol );
        Gemv( NORMAL, Field(1), pinvA, b, q );
    }
    else if( m >= n )
    {
        Q = A;
        qr::Explicit( Q, R );
        Gemv( ADJOINT, Field(1), Q, b, q );
        Trsv( UPPER, NORMAL, NON_UNIT, R, q );
    }
    else
    {
        Q = A;
        lq::Explicit( L, Q );
        s = b;
        Trsv( LOWER, NORMAL, NON_UNIT, L, s );
        Gemv( ADJOINT, Field(1), Q, s, q );
    }

    if( ctrl.progress )
    {
        const Real qOneNorm = OneNorm( q );
        Output(" || pinv(A) b ||_1 = ",qOneNorm);
    }

    // Start the basis pursuit
    Int numIter=0;
    Matrix<Field> x, u, t, zOld, xHat;
    Zeros( x, n, 1 );
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    while( numIter < ctrl.maxIter )
    {
        zOld = z;

        // x := P*(z-u) + q
        //    = (I-pinv(A)*A)(z-u) + q
        //    = (z-u) - pinv(A)*A*(z-u) + q
        s = z;
        s -= u;
        x = s;
        Gemv( NORMAL, Field(1), A, s, t );
        if( ctrl.usePinv )
        {
            Gemv( NORMAL, Field(1), pinvA, t, s );
        }
        else if( m >= n )
        {
            Gemv( ADJOINT, Field(1), Q, t, s );
            Trsv( UPPER, NORMAL, NON_UNIT, R, s );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, L, t );
            Gemv( ADJOINT, Field(1), Q, t, s );
        }
        x -= s;
        x += q;

        // xHat := alpha x + (1-alpha) zOld
        xHat = x;
        xHat *= ctrl.alpha;
        Axpy( 1-ctrl.alpha, zOld, xHat );

        // z := SoftThresh(xHat+u,1/rho)
        z = xHat;
        z += u;
        SoftThreshold( z, 1/ctrl.rho );

        // u := u + (xHat - z)
        u += xHat;
        u -= z;

        // rNorm := || x - z ||_2
        s = x;
        s -= z;
        const Real rNorm = FrobeniusNorm( s );

        // sNorm := || rho*(z-zOld) ||_2
        s = z;
        s -= zOld;
        const Real sNorm = Abs(ctrl.rho)*FrobeniusNorm( s );

        const Real epsPri = Sqrt(Real(n))*ctrl.absTol +
            ctrl.relTol*Max(FrobeniusNorm(x),FrobeniusNorm(z));
        const Real epsDual = Sqrt(Real(n))*ctrl.absTol +
            ctrl.relTol*Abs(ctrl.rho)*FrobeniusNorm(u);

        if( ctrl.progress )
        {
            const Real xOneNorm = OneNorm( x );
            Output
            (numIter,": ||x-z||_2=",rNorm,", epsPri=",epsPri,
             ", |rho| ||z-zOld||_2=",sNorm,", and epsDual=",epsDual,
             ", ||x||_1=",xOneNorm);
        }

        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( ctrl.maxIter == numIter )
        RuntimeError("Basis pursuit failed to converge");
    return numIter;
}

template<typename Field>
Int ADMM
( const AbstractDistMatrix<Field>& APre,
  const AbstractDistMatrix<Field>& bPre,
        AbstractDistMatrix<Field>& zPre,
  const ADMMCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Field,Field,MC,MR>
      AProx( APre ),
      bProx( bPre );
    DistMatrixWriteProxy<Field,Field,MC,MR>
      zProx( zPre );
    auto& A = AProx.GetLocked();
    auto& b = bProx.GetLocked();
    auto& z = zProx.Get();

    // Find a means of quickly applyinv pinv(A) and then form pinv(A) b
    // NOTE: If m >= n and A has full column rank, then basis pursuit is
    //       irrelevant, as there is a unique solution, which is found
    //       through least squares. If A does *not* have full column rank,
    //       then the QR factorization is not enough.
    //       For the same reason, the LQ factorization will fail if A does
    //       not have full row rank.
    //
    // TODO(poulson): Instead form a basis for the null-space of A
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    DistMatrix<Field> pinvA(grid);
    DistMatrix<Field> Q(grid), L(grid), R(grid);
    DistMatrix<Field> q(grid), s(grid);
    if( ctrl.usePinv )
    {
        pinvA = A;
        Pseudoinverse( pinvA, ctrl.pinvTol );
        Gemv( NORMAL, Field(1), pinvA, b, q );
    }
    else if( m >= n )
    {
        Q = A;
        qr::Explicit( Q, R );
        Gemv( ADJOINT, Field(1), Q, b, q );
        Trsv( UPPER, NORMAL, NON_UNIT, R, q );
    }
    else
    {
        Q = A;
        lq::Explicit( L, Q );
        s = b;
        Trsv( LOWER, NORMAL, NON_UNIT, L, s );
        Gemv( ADJOINT, Field(1), Q, s, q );
    }

    if( ctrl.progress )
    {
        const Real qOneNorm = OneNorm( q );
        if( grid.Rank() == 0 )
            Output(" || pinv(A) b ||_1 = ",qOneNorm);
    }

    // Start the basis pursuit
    Int numIter=0;
    DistMatrix<Field> x(grid), u(grid), t(grid), zOld(grid), xHat(grid);
    Zeros( x, n, 1 );
    Zeros( z, n, 1 );
    Zeros( u, n, 1 );
    while( numIter < ctrl.maxIter )
    {
        zOld = z;

        // x := P*(z-u) + q
        //    = (I-pinv(A)*A)(z-u) + q
        //    = (z-u) - pinv(A)*A*(z-u) + q
        s = z;
        s -= u;
        x = s;
        Gemv( NORMAL, Field(1), A, s, t );
        if( ctrl.usePinv )
        {
            Gemv( NORMAL, Field(1), pinvA, t, s );
        }
        else if( m >= n )
        {
            Gemv( ADJOINT, Field(1), Q, t, s );
            Trsv( UPPER, NORMAL, NON_UNIT, R, s );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, L, t );
            Gemv( ADJOINT, Field(1), Q, t, s );
        }
        x -= s;
        x += q;

        // xHat := alpha x + (1-alpha) zOld
        xHat = x;
        xHat *= ctrl.alpha;
        Axpy( 1-ctrl.alpha, zOld, xHat );

        // z := SoftThresh(xHat+u,1/rho)
        z = xHat;
        z += u;
        SoftThreshold( z, 1/ctrl.rho );

        // u := u + (xHat - z)
        u += xHat;
        u -= z;

        // rNorm := || x - z ||_2
        s = x;
        s -= z;
        const Real rNorm = FrobeniusNorm( s );

        // sNorm := || rho*(z-zOld) ||_2
        s = z;
        s -= zOld;
        const Real sNorm = Abs(ctrl.rho)*FrobeniusNorm( s );

        const Real epsPri = Sqrt(Real(n))*ctrl.absTol +
            ctrl.relTol*Max(FrobeniusNorm(x),FrobeniusNorm(z));
        const Real epsDual = Sqrt(Real(n))*ctrl.absTol +
            ctrl.relTol*Abs(ctrl.rho)*FrobeniusNorm(u);

        if( ctrl.progress )
        {
            const Real xOneNorm = OneNorm( x );
            if( grid.Rank() == 0 )
                Output
                (numIter,": ||x-z||_2=",rNorm,", epsPri=",epsPri,
                 ", |rho| ||z-zOld||_2=",sNorm,", and epsDual=",epsDual,
                 ", ||x||_1=",xOneNorm);
        }

        if( rNorm < epsPri && sNorm < epsDual )
            break;
        ++numIter;
    }
    if( ctrl.maxIter == numIter )
        RuntimeError("Basis pursuit failed to converge");
    return numIter;
}

} // namespace bp
} // namespace El
