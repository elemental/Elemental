/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;
    const int commRank = El::mpi::Rank( comm );
    const int commSize = El::mpi::Size( comm );

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int n = El::Input("--size","size of matrix",100);
        const El::Int numRhs = El::Input("--numRhs","# of right-hand sides",1);
        const El::Int blocksize =
          El::Input("--blocksize","algorithmic blocksize",64);
        const El::Int numTests = El::Input("--numTests","number of tests",3);
        const bool error = El::Input("--error","test Elemental error?",true);
        El::Int gridHeight = El::Input("--gridHeight","grid height",0);
        const bool details = El::Input("--details","print norm details?",false);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::SetBlocksize( blocksize );

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = El::Grid::DefaultHeight( commSize );
        El::Grid grid( comm, gridHeight );
        if( commRank == 0 )
            El::Output("Grid is: ",grid.Height()," x ",grid.Width());

        // Set up random A and B, then make the copies X := B
        El::Timer timer;
        El::DistMatrix<Scalar> A(grid), B(grid), X(grid);

        for( El::Int test=0; test<numTests; ++test )
        {
            El::Uniform( A, n, n );
            El::Uniform( B, n, numRhs );
            X = B;
            if( print )
            {
                El::Print( A, "A" );
                El::Print( B, "B" );
            }

            // Perform the LU factorization and simultaneous solve
            if( commRank == 0 )
                El::Output("Starting Elemental linear solve");
            El::mpi::Barrier( comm );
            if( commRank == 0 )
                timer.Start();
            El::LinearSolve( A, X );
            El::mpi::Barrier( comm );
            if( commRank == 0 )
                El::Output(timer.Stop()," seconds");

            if( error )
            {
                // Form R := A X - B
                auto R( B );
                El::Gemm
                ( El::NORMAL, El::NORMAL, Scalar(1), A, X, Scalar(-1), R );

                // Compute infinity norms and a relative residual
                const Real eps = El::limits::Epsilon<Real>();
                const Real AInfNorm = El::InfinityNorm( A );
                const Real BInfNorm = El::InfinityNorm( B );
                const Real XInfNorm = El::InfinityNorm( X );
                const Real RInfNorm = El::InfinityNorm( R );
                const Real infResidual = RInfNorm / (AInfNorm*XInfNorm*eps*n);
                if( commRank == 0 )
                {
                    if( details )
                    {
                        El::Output("");
                        El::Output("||A||_oo       = ",AInfNorm);
                        El::Output("||B||_oo       = ",BInfNorm);
                        El::Output("||X||_oo       = ",XInfNorm);
                        El::Output("||A X - B||_oo = ",RInfNorm);
                    }
                    El::Output
                    ("||A X - B||_oo / (||A||_oo ||X||_oo eps n) = ",
                     infResidual);
                }

                // Compute one norms and a relative residual
                const Real AOneNorm = El::OneNorm( A );
                const Real BOneNorm = El::OneNorm( B );
                const Real XOneNorm = El::OneNorm( X );
                const Real ROneNorm = El::OneNorm( R );
                const Real oneResidual = ROneNorm / (AOneNorm*XOneNorm*eps*n);
                if( commRank == 0 )
                {
                    if( details )
                    {
                        El::Output("");
                        El::Output("||A||_1       = ",AOneNorm);
                        El::Output("||B||_1       = ",BOneNorm);
                        El::Output("||X||_1       = ",XOneNorm);
                        El::Output("||A X - B||_1 = ",ROneNorm);
                    }
                    El::Output
                    ("||A X - B||_1 / (||A||_1 ||X||_1 eps n) = ",
                     oneResidual,"\n");
                }
            }
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
