/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_UNIFORM_INC
using namespace El;

// This driver calls an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/lasso/lasso.html
//
// The Least Absolute Shrinkage and Selection Operator (LASSO)
//   minimizes || A x - b ||_2^2 + lambda || x ||_1,
// which is equivalent to minimizing || A x - b ||_2 subject to 
// || x ||_1 <= t for some t >= 0.

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
        const Real lambda = Input("--lambda","one-norm coefficient",1.);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = Input("--relTol","relative tolerance",1e-4);
        const bool inv = Input("--inv","use explicit inverse",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A, b;
        Uniform( A, m, n );
        Uniform( b, m, 1 );
        if( print )
        {
            Print( A,     "A"     );
            Print( b,     "b"     );
        }
        if( display )
            Display( A, "A" );

        DistMatrix<C> z;
        Lasso
        ( A, b, lambda, z, rho, alpha, maxIter, absTol, relTol, inv, progress );
        if( print )
            Print( z, "z" );
        const Real zOneNorm = OneNorm( z );
        const Real zTwoNorm = FrobeniusNorm( z );
        const Int  zZeroNorm = ZeroNorm( z );
        const Real bTwoNorm = FrobeniusNorm( b );
        Gemv( NORMAL, C(-1), A, z, C(1), b );
        const Real rTwoNorm = FrobeniusNorm( b );
        if( mpi::Rank(mpi::COMM_WORLD) == 0 )
        {
            std::cout << "|| A z - b ||_2 = " << rTwoNorm << "\n"
                      << "|| b ||_2 = " << bTwoNorm << "\n"
                      << "|| z ||_1 = " << zOneNorm << "\n"
                      << "|| z ||_0 = " << zZeroNorm << "\n" << std::endl;
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
