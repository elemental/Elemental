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

// Solve
//
//     minimize || A z - y ||_2 such that z >= 0
//        z 
//
// via the Quadratic Program
//
//     minimize    (1/2) x' P x + q' x 
//     subject to  x >= 0
//
// with P = A^T A and q = -A^H y.
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

        DistMatrix<Real> A, y;
        Uniform( A, m, n );
        Uniform( y, m, 1 );
        if( print )
        {
            Print( A, "A" );
            Print( y, "y" );
        }
        if( display )
            Display( A, "A" );

        DistMatrix<Real> z;
        NonNegativeLeastSquares
        ( A, y, z, rho, alpha, maxIter, absTol, relTol, inv, progress );

        if( print )
            Print( z, "z" );

        const double yNorm = FrobeniusNorm( y );
        Gemv( NORMAL, Real(-1), A, z, Real(1), y );
        const double eNorm = FrobeniusNorm( y );
        if( mpi::WorldRank() == 0 )
            std::cout << "|| y - A z ||_2 / || y ||_2 = " << eNorm/yNorm
                      << std::endl;
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
