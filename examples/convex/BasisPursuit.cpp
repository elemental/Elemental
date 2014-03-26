/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_GEMV_INC
#include ELEM_TRSV_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_LQ_INC
#include ELEM_QR_INC
#include ELEM_PSEUDOINVERSE_INC
#include ELEM_SOFTTHRESHOLD_INC
#include ELEM_UNIFORM_INC
#include ELEM_ZEROS_INC
using namespace elem;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/basis_pursuit/basis_pursuit.html
// which is derived from the distributed ADMM article of Boyd et al.

typedef double Real;
typedef Complex<Real> C;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int m = Input("--m","height of matrix",100);
        const Int n = Input("--n","width of matrix",200);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real probNnz = Input("--probNnz","prob. of nonzero x entry",0.1);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1.e-4);
        const Real relTol = Input("--relTol","relative tolerance",1.e-2);
        const Real pinvTol = Input("--pinvTol","pseudoinv tolerance",0.);
        const bool usePinv = Input("--usePinv","Directly compute pinv(A)",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A, b, xTrue;
        Uniform( A, m, n );
        Zeros( xTrue, n, 1 );
        if( xTrue.LocalWidth() == 1 )
        {
            for( Int iLoc=0; iLoc<xTrue.LocalHeight(); ++iLoc )
                if( SampleUniform<Real>() <= probNnz )
                    xTrue.SetLocal( iLoc, 0, SampleBall<C>() );
        }
        Gemv( NORMAL, C(1), A, xTrue, b ); 
        if( print )
        {
            Print( A,     "A"     );
            Print( xTrue, "xTrue" );
            Print( b,     "b"     );
        }
        if( display )
            Display( A, "A" );

        // Find a means of quickly applyinv pinv(A) and then form pinv(A) b
        // NOTE: If m >= n and A has full column rank, then basis pursuit is 
        //       irrelevant, as there is a unique solution, which is found 
        //       through least squares. If A does *not* have full column rank,
        //       then the QR factorization is not enough.
        //       For the same reason, the LQ factorization will fail if A does
        //       not have full row rank.
        DistMatrix<C> pinvA;
        DistMatrix<C> Q, L, R;
        DistMatrix<C> q, s;
        if( usePinv )
        {
            pinvA = A;
            Pseudoinverse( pinvA, pinvTol );
            Gemv( NORMAL, C(1), pinvA, b, q );
        }
        else if( m >= n )
        {
            Q = A; 
            qr::Explicit( Q, R );
            Gemv( ADJOINT, C(1), Q, b, q );
            Trsv( UPPER, NORMAL, NON_UNIT, R, q );
        }
        else
        {
            Q = A;
            lq::Explicit( L, Q );
            s = b;
            Trsv( LOWER, NORMAL, NON_UNIT, L, s );
            Gemv( ADJOINT, C(1), Q, s, q );
        }
        if( print )
            Print( q, "q := pinv(A) b" );

        const Real xTrueOneNorm = OneNorm( xTrue );
        const Real qOneNorm = OneNorm( q );
        if( mpi::WorldRank() == 0 )
        {
            std::cout << " || xTrue     ||_1 = " << xTrueOneNorm << "\n"
                      << " || pinv(A) b ||_1 = " << qOneNorm << std::endl;
        }

        // Start the basis pursuit
        Int numIter=0;
        DistMatrix<C> x, z, u, t;
        DistMatrix<C> zOld, xHat;
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
            Axpy( C(-1), u, s );
            x = s;
            Gemv( NORMAL, C(1), A, s, t );
            if( usePinv )
            {
                Gemv( NORMAL, C(1), pinvA, t, s );
            }
            else if( m >= n )
            {
                Gemv( ADJOINT, C(1), Q, t, s );
                Trsv( UPPER, NORMAL, NON_UNIT, R, s );
            }
            else
            {
                Trsv( LOWER, NORMAL, NON_UNIT, L, t );
                Gemv( ADJOINT, C(1), Q, t, s );
            }
            Axpy( C(-1), s, x );
            Axpy( C(1),  q, x );

            // xHat := alpha x + (1-alpha) zOld
            xHat = x;
            Scale( alpha, xHat );
            Axpy( 1-alpha, zOld, xHat );

            // z := SoftThresh(xHat+u,1/rho)
            z = xHat;
            Axpy( C(1), u, z );
            SoftThreshold( z, 1/rho );

            // u := u + (xHat - z)
            Axpy( C(1),  xHat, u );
            Axpy( C(-1), z,    u );

            // rNorm := || x - z ||_2
            s = x;
            Axpy( C(-1), z, s );
            const Real rNorm = FrobeniusNorm( s );

            // sNorm := || rho*(z-zOld) ||_2
            s = z;
            Axpy( C(-1), zOld, s );
            const Real sNorm = Abs(rho)*FrobeniusNorm( s );
            
            const Real epsPri = Sqrt(Real(n))*absTol + 
                relTol*Max(FrobeniusNorm(x),FrobeniusNorm(z));
            const Real epsDual = Sqrt(Real(n))*absTol +
                relTol*Abs(rho)*FrobeniusNorm(u);

            if( progress )
            {
                const Real xOneNorm = OneNorm( x );
                if( mpi::WorldRank() == 0 )
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
        if( maxIter == numIter && mpi::WorldRank() == 0 )
            std::cout << "Basis pursuit failed to converge" << std::endl;
        if( print )
        {
            Print( x, "x" );
            Print( z, "z" );
            Print( u, "u" );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
