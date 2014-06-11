/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_QUADRATICPROGRAM_INC
#include EL_HERMITIANUNIFORMSPECTRUM_INC
#include EL_GAUSSIAN_INC
using namespace El;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This example attempts to solve the following convex optimization problem:
//     minimize    (1/2) x' P x + q' x 
//     subject to  lb <= x <= ub
//

typedef double Real;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int m = Input("--m","matrix height",200);
        const Int n = Input("--n","problem size",100);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real lb = Input("--lb","lower bound for x",0.5);
        const Real ub = Input("--ub","upper bound for x",1.0);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = Input("--relTol","relative tolerance",1e-4);
        const bool inv = Input("--inv","form inv(LU) to avoid trsv?",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Real> A, P, q, y;
        Uniform( A, m, n );
        Herk( LOWER, ADJOINT, 1., A, P );
        Uniform( y, m, 1 );
        Gemv( ADJOINT, 1., A, y, q );
        if( print )
        {
            Print( A, "A" );
            Print( P, "P" );
            Print( y, "y" );
            Print( q, "q" );
        }
        if( display )
        {
            Display( A, "A" );
            Display( P, "P" );
        }

        DistMatrix<Real> x, z, u;
        QuadraticProgram
        ( P, q, 0., 1e6, x, z, u, rho, alpha, maxIter, absTol, relTol, inv, 
          progress );

        if( print )
        {
            Print( x, "x" );
            Print( z, "z" );
            Print( u, "u" );
        }

        const double yNorm = FrobeniusNorm( y );
        Gemv( NORMAL, -1., A, x, 1., y );
        const double eNorm = FrobeniusNorm( y );
        if( mpi::WorldRank() == 0 )
            std::cout << "|| y - A x ||_2 / || y ||_2 = " << eNorm/yNorm
                      << std::endl;
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
