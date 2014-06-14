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

// This driver uses an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/linprog/linprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// It attempts to solve the following linear program:
//     minimize    c' x
//     subject to  A x = b, x >= 0
//

typedef double Real;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int m = Input("--m","height of matrix",100);
        const Int n = Input("--n","width of matrix",200);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-4);
        const Real relTol = Input("--relTol","relative tolerance",1e-2);
        const bool inv = Input("--inv","form inv(LU) to avoid trsv?",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Real> A, b, c, xTrue;
        Uniform( A, m, n, 1., 1. ); // mean=radius=1, so sample in [0,2]
        Zeros( xTrue, n, 1 );
        if( xTrue.LocalWidth() == 1 )
            for( Int iLoc=0; iLoc<xTrue.LocalHeight(); ++iLoc )
                xTrue.SetLocal( iLoc, 0, Abs(SampleNormal<Real>()) );
        Gemv( NORMAL, Real(1), A, xTrue, b ); 
        Uniform( c, n, 1, 1., 1. ); // mean=radius=1, so sample in [0,2]
        if( print )
        {
            Print( A,     "A"     );
            Print( xTrue, "xTrue" );
            Print( b,     "b"     );
            Print( c,     "c"     );
        }
        if( display )
            Display( A, "A" );
        const Real objectiveTrue = Dot( c, xTrue );
        if( mpi::WorldRank() == 0 )
            std::cout << "c'xTrue=" << objectiveTrue << std::endl;

        DistMatrix<Real> z;
        LinearProgram
        ( A, b, c, z, rho, alpha, maxIter, absTol, relTol, inv, progress );

        if( print )
            Print( z, "z" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
