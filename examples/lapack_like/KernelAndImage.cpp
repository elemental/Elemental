/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int rank = Input("--rank","rank of matrix",10);
        const Int blocksize = Input("--blocksize","algorithmic blocksize",32);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( blocksize );

        const int commRank = mpi::Rank();
        Timer timer;

        Grid g( mpi::COMM_WORLD );
        if( commRank == 0 )
            Output("Grid is ",g.Height()," x ",g.Width());
        DistMatrix<C> A(g), X(g), Y(g);
        Uniform( X, m, rank );
        Uniform( Y, rank, n );
        Gemm( NORMAL, NORMAL, C(1), X, Y, A );
        if( print )
            Print( A, "A" );

        // Compute the image and kernel
        if( commRank == 0 )
            timer.Start();
        DistMatrix<C> M(g), K(g);
        ImageAndKernel( A, M, K );
        if( commRank == 0 )
            Output("ImageAndKernel took ",timer.Stop()," seconds");
        if( print )
        {
            Print( M, "M" );
            Print( K, "K" );
        }

        // Compute just the image
        if( commRank == 0 )
            timer.Start();
        DistMatrix<C> MOnly(g);
        Image( A, MOnly );
        if( commRank == 0 )
            Output("Image took ",timer.Stop()," seconds");
        if( print )
            Print( MOnly, "MOnly" );

        // Compute just the kernel
        if( commRank == 0 )
            timer.Start();
        DistMatrix<C> KOnly(g);
        Kernel( A, KOnly );
        if( commRank == 0 )
            Output("Kernel took ",timer.Stop()," seconds");
        if( print )
            Print( KOnly, "KOnly" );

        // TODO: Add correctness tests 
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
