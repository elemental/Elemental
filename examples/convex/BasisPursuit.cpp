/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_ZERONORM_INC
#include ELEM_BASISPURSUIT_INC
#include ELEM_UNIFORM_INC
using namespace elem;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/basis_pursuit/basis_pursuit.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// Basis pursuit seeks the solution to A x = b which minimizes || x ||_1

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
        const Real absTol = Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = Input("--relTol","relative tolerance",1e-4);
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
        const Int trueZeroNorm = ZeroNorm( xTrue );
        Gemv( NORMAL, C(1), A, xTrue, b ); 
        if( print )
        {
            Print( A,     "A"     );
            Print( xTrue, "xTrue" );
            Print( b,     "b"     );
        }
        if( display )
            Display( A, "A" );

        DistMatrix<C> x, z, u;
        BasisPursuit
        ( A, b, x, z, u, alpha, rho, maxIter, absTol, relTol, usePinv, pinvTol,
          progress );
        if( print )
        {
            Print( x, "x" );
            Print( z, "z" );
            Print( u, "u" );
        }
        const Int xZeroNorm = ZeroNorm( x );
        const Int zZeroNorm = ZeroNorm( z );
        const Int uZeroNorm = ZeroNorm( u );
        if( mpi::Rank(mpi::COMM_WORLD) == 0 )
        {
            std::cout << "|| xTrue ||_0 = " << trueZeroNorm << "\n"
                      << "|| x     ||_0 = " << xZeroNorm << "\n"
                      << "|| z     ||_0 = " << zZeroNorm << "\n" 
                      << "|| u     ||_0 = " << uZeroNorm << "\n" << std::endl;
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
