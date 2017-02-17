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

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int m = El::Input("--height","height of matrix",100);
        const El::Int n = El::Input("--width","width of matrix",100);
        const El::Int rank = El::Input("--rank","rank of matrix",10);
        const El::Int blocksize =
          El::Input("--blocksize","algorithmic blocksize",32);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::SetBlocksize( blocksize );

        El::mpi::Comm comm = El::mpi::COMM_WORLD;
        const int commRank = El::mpi::Rank(comm);
        El::Timer timer;

        El::Grid grid( comm );
        if( commRank == 0 )
            El::Output("Grid is ",grid.Height()," x ",grid.Width());
        El::DistMatrix<Scalar> A(grid), X(grid), Y(grid);
        El::Uniform( X, m, rank );
        El::Uniform( Y, rank, n );
        El::Gemm( El::NORMAL, El::NORMAL, Scalar(1), X, Y, A );
        if( print )
            El::Print( A, "A" );

        // Compute the image and kernel
        if( commRank == 0 )
            timer.Start();
        El::DistMatrix<Scalar> M(grid), K(grid);
        El::ImageAndKernel( A, M, K );
        if( commRank == 0 )
            El::Output("ImageAndKernel took ",timer.Stop()," seconds");
        if( print )
        {
            El::Print( M, "M" );
            El::Print( K, "K" );
        }

        // Compute just the image
        if( commRank == 0 )
            timer.Start();
        El::DistMatrix<Scalar> MOnly(grid);
        El::Image( A, MOnly );
        if( commRank == 0 )
            El::Output("Image took ",timer.Stop()," seconds");
        if( print )
            El::Print( MOnly, "MOnly" );

        // Compute just the kernel
        if( commRank == 0 )
            timer.Start();
        El::DistMatrix<Scalar> KOnly(grid);
        El::Kernel( A, KOnly );
        if( commRank == 0 )
            El::Output("Kernel took ",timer.Stop()," seconds");
        if( print )
            El::Print( KOnly, "KOnly" );

        // TODO(poulson): Add correctness tests
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
