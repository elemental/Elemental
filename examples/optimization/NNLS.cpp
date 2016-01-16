/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// Solve
//
//     minimize || A z - y ||_2 such that z >= 0
//        z 
//
// via the Quadratic Program
//
//     minimize    (1/2) x' Q x + c' x 
//     subject to  x >= 0
//
// with Q = A^T A and c = -A^H y.
//

typedef double Real;

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int m = Input("--m","matrix height",150);
        const Int n = Input("--n","matrix width",100);
        const Int k = Input("--k","number of right-hand sides",2);
        // TODO: Test both the ADMM and IPM versions
        /*
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = Input("--relTol","relative tolerance",1e-4);
        const bool inv = Input("--inv","form inv(LU) to avoid trsv?",true);
        const bool progress = Input("--progress","print progress?",true);
        */
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        const Real center = 1;
        const Real radius = 1;

        DistMatrix<Real> A, B;
        Uniform( A, m, n, center, radius );
        Uniform( B, m, k, center, radius );
        if( print )
            Print( A, "A" );
        if( display )
            Display( A, "A" );

        /*
        qp::box::ADMMCtrl<Real> ctrl;
        ctrl.rho = rho;
        ctrl.alpha = alpha;
        ctrl.maxIter = maxIter;
        ctrl.absTol = absTol;
        ctrl.relTol = relTol;
        ctrl.inv = inv;
        ctrl.print = progress;

        DistMatrix<Real> X;
        nnls::ADMM( A, B, X, ctrl );
        if( print )
            Print( X, "X" );
        */

        Timer timer;
        DistMatrix<Real> X;
        if( mpi::Rank() == 0 )
            timer.Start();
        NNLS( A, B, X );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( print )
        {
            Print( B, "B" );
            Print( X, "X" );
        }
        if( display )
        {
            Display( B, "B" );
            Display( X, "X" );
        }

        const Real ANorm = FrobeniusNorm( A );
        Gemm( NORMAL, TRANSPOSE, Real(-1), B, X, Real(1), A );
        const Real ENorm = FrobeniusNorm( A );
        if( mpi::Rank() == 0 )
        {
            Output("NNLS time: ",timer.Total(),"secs");
            Output("|| A - B X' ||_2 / || A ||_2 = ",ENorm/ANorm);
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
