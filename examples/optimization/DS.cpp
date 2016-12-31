/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

typedef double Real;

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int m = El::Input("--m","height of matrix",100);
        const El::Int n = El::Input("--n","width of matrix",200);
        // TODO(poulson): Add options for controlling IPM
        const El::Int maxIter =
          El::Input("--maxIter","maximum # of iter's",500);
        const Real lambda = El::Input("--lambda","DS parameter",0.5);
        const bool display = El::Input("--display","display matrices?",false);
        const bool print = El::Input("--print","print matrices",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::SparseMatrix<Real> A;
        El::Matrix<Real> b, xTrue;
        El::Identity( A, m, n );
        El::Uniform( b, m, 1 );
        if( print )
        {
            El::Print( A, "A" );
            El::Print( b, "b" );
        }
        if( display )
            El::Display( A, "A" );

        El::lp::affine::Ctrl<Real> affineCtrl; 
        affineCtrl.mehrotraCtrl.print = true;

        El::Matrix<Real> x;
        El::Timer timer;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::DS( A, b, lambda, x, affineCtrl );
        if( El::mpi::Rank() == 0 )
            timer.Stop();
        if( print )
            El::Print( x, "x" );
        const El::Int xZeroNorm = El::ZeroNorm( x );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("Dantzig Selector time: ",timer.Total()," secs");
            El::Output("|| x ||_0 = ",xZeroNorm);
        }
        SoftThreshold( x, El::Sqrt(El::limits::Epsilon<Real>()) );
        if( print )
            El::Print( x, "xThresh" );
        const El::Int xZeroNormThresh = El::ZeroNorm( x );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("|| xThresh ||_0 = ",xZeroNormThresh);
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
