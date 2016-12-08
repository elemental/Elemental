/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int m = El::Input("--height","height of matrix",100);
        const El::Int n = El::Input("--width","width of matrix",100);
        const El::SignScaling scaling =
          static_cast<El::SignScaling>
          (El::Input("--scaling","scaling strategy",0));
        const El::Int maxIts = El::Input("--maxIts","max number of iter's",100);
        const double tol = El::Input("--tol","convergence tolerance",1e-6);
        const bool progress =
          El::Input("--progress","print sign progress?",true);
        const bool print = El::Input("--print","print matrix?",false);
        const bool display = El::Input("--display","display matrix?",false);
        El::ProcessInput();
        El::PrintInputReport();

        const El::Grid grid( comm );
        El::DistMatrix<Scalar> A(grid);
        El::Uniform( A, m, n );
        if( print )
            El::Print( A, "A" );
        if( display )
            El::Display( A, "A" );

        El::SignCtrl<Real> signCtrl;
        signCtrl.maxIts = maxIts;
        signCtrl.tol = tol;
        signCtrl.progress = progress;
        signCtrl.scaling = scaling;

        El::Timer timer;
        // Compute sgn(A)
        if( El::mpi::Rank(comm) == 0 )
            timer.Start();
        El::Sign( A, signCtrl );
        if( El::mpi::Rank(comm) == 0 )
            timer.Stop();
        if( print )
            El::Print( A, "A" );
        if( display )
            El::Display( A, "A" );
        if( El::mpi::Rank(comm) == 0 )
            El::Output("Sign time: ",timer.Total()," secs");
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
