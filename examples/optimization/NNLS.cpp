/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

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
    El::Environment env( argc, argv );

    try
    {
        const El::Int m = El::Input("--m","matrix height",150);
        const El::Int n = El::Input("--n","matrix width",100);
        const El::Int k = El::Input("--k","number of right-hand sides",2);
        // TODO(poulson): Test both the ADMM and IPM versions
        /*
        const El::Int maxIter =
          El::Input("--maxIter","maximum # of iter's",500);
        const Real rho = El::Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = El::Input("--alpha","over-relaxation",1.2);
        const Real absTol = El::Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = El::Input("--relTol","relative tolerance",1e-4);
        const bool inv = El::Input("--inv","form inv(LU) to avoid trsv?",true);
        const bool progress = El::Input("--progress","print progress?",true);
        */
        const bool display = El::Input("--display","display matrices?",false);
        const bool print = El::Input("--print","print matrices",false);
        El::ProcessInput();
        El::PrintInputReport();

        const Real center = 1;
        const Real radius = 1;

        El::DistMatrix<Real> A, B;
        El::Uniform( A, m, n, center, radius );
        El::Uniform( B, m, k, center, radius );
        if( print )
            El::Print( A, "A" );
        if( display )
            El::Display( A, "A" );

        /*
        El::qp::box::ADMMCtrl<Real> ctrl;
        ctrl.rho = rho;
        ctrl.alpha = alpha;
        ctrl.maxIter = maxIter;
        ctrl.absTol = absTol;
        ctrl.relTol = relTol;
        ctrl.inv = inv;
        ctrl.print = progress;

        El::DistMatrix<Real> X;
        El::nnls::ADMM( A, B, X, ctrl );
        if( print )
            El::Print( X, "X" );
        */

        El::Timer timer;
        El::DistMatrix<Real> X;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::NNLS( A, B, X );
        if( El::mpi::Rank() == 0 )
            timer.Stop();
        if( print )
        {
            El::Print( B, "B" );
            El::Print( X, "X" );
        }
        if( display )
        {
            El::Display( B, "B" );
            El::Display( X, "X" );
        }

        const Real ANorm = El::FrobeniusNorm( A );
        El::Gemm( El::NORMAL, El::TRANSPOSE, Real(-1), B, X, Real(1), A );
        const Real ENorm = El::FrobeniusNorm( A );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("NNLS time: ",timer.Total()," secs");
            El::Output("|| A - B X' ||_2 / || A ||_2 = ",ENorm/ANorm);
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
