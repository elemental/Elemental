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
    const int commSize = El::mpi::Size( comm );

    try
    {
        const El::Int m = El::Input("--m","height of A",100);
        const El::Int n = El::Input("--n","width of A",100);
        const El::Int p = El::Input("--p","width of B",20);
        const El::Int numRhs = El::Input("--numRhs","# of right-hand sides",5);
        const El::Int blocksize =
          El::Input("--blocksize","algorithmic blocksize",64);
        const bool print = El::Input("--print","print matrices?",false);
        El::Int gridHeight = El::Input("--gridHeight","grid height",0);
        El::ProcessInput();
        El::PrintInputReport();

        // Set the algorithmic blocksize
        El::SetBlocksize( blocksize );

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = El::Grid::DefaultHeight( commSize );
        El::Grid grid( comm, gridHeight );

        El::DistMatrix<El::Complex<double>>
          A(grid), B(grid), D(grid), X(grid), Y(grid);
        El::Uniform( A, m, n );
        El::Uniform( B, m, p );
        El::Uniform( D, m, numRhs );
        if( print )
        {
            El::Print( A, "A" );
            El::Print( B, "B" );
            El::Print( D, "D" );
        }

        El::Timer timer;
        if( El::mpi::Rank(comm) == 0 )
            timer.Start();
        El::GLM( A, B, D, X, Y );
        if( El::mpi::Rank(comm) == 0 )
            timer.Stop();
        if( print )
        {
            El::Print( X, "X" );
            El::Print( Y, "Y" );
        }

        const double DFrob = El::FrobeniusNorm( D );
        El::Gemm
        ( El::NORMAL, El::NORMAL,
          El::Complex<double>(-1), A, D, El::Complex<double>(1), D );
        El::Gemm
        ( El::NORMAL, El::NORMAL,
          El::Complex<double>(-1), B, Y, El::Complex<double>(1), D );
        const double EFrob = El::FrobeniusNorm( D );
        if( print )
            El::Print( D, "D - A X - B Y" );
        if( El::mpi::Rank(comm) == 0 )
        {
            El::Output("GLM time: ",timer.Total()," secs");
            El::Output
            ("|| D             ||_F = ",DFrob,"\n",
             "|| A X + B Y - D ||_F = ",EFrob,"\n");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
