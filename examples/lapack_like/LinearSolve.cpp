/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try 
    {
        const Int n = Input("--size","size of matrix",100);
        const Int numRhs = Input("--numRhs","# of right-hand sides",1); 
        const Int blocksize = Input("--blocksize","algorithmic blocksize",64);
        const Int numTests = Input("--numTests","number of tests",3);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true); 
#else
        const bool scalapack = false;
#endif
        const Int mb = Input("--mb","block height",32);
        const Int nb = Input("--nb","block width",32);
        Int gridHeight = Input("--gridHeight","grid height",0);
        const bool details = Input("--details","print norm details?",false);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( blocksize );
        SetDefaultBlockHeight( mb ); 
        SetDefaultBlockWidth( nb );

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( commSize );
        Grid grid( comm, gridHeight );
        if( commRank == 0 )
            Output("Grid is: ",grid.Height()," x ",grid.Width());

        // Set up random A and B, then make the copies X := B
        Timer timer;
        DistMatrix<double> A(grid), B(grid), X(grid);
        for( Int test=0; test<numTests; ++test )
        {
            Uniform( A, n, n );
            Uniform( B, n, numRhs );
            X = B;
            if( print )
            {
                Print( A, "A" );
                Print( B, "B" );
            }

            if( scalapack )
            {
                if( commRank == 0 )
                {
                    cout << "Starting ScaLAPACK linear solve...";
                    cout.flush();
                }
                DistMatrix<double,MC,MR,BLOCK> ABlock( A ), BBlock( B );
                mpi::Barrier( comm );
                if( commRank == 0 )
                    timer.Start();
                LinearSolve( ABlock, BBlock );
                if( commRank == 0 )
                    Output(timer.Stop()," seconds");
            }

            // Perform the LU factorization and simultaneous solve
            if( commRank == 0 )
            {
                cout << "Starting linear solve...";
                cout.flush();
            }
            mpi::Barrier( comm );
            if( commRank == 0 )
                timer.Start();
            LinearSolve( A, X );
            mpi::Barrier( comm );
            if( commRank == 0 )
                Output(timer.Stop()," seconds");

            // Form R := A X - B
            DistMatrix<> R( B );
            Gemm( NORMAL, NORMAL, 1., A, X, -1., R );

            // Compute the relevant infinity norms and a relative residual
            const double epsilon = lapack::MachineEpsilon();
            const double AInfNorm = InfinityNorm( A );
            const double BInfNorm = InfinityNorm( B );
            const double XInfNorm = InfinityNorm( X );
            const double RInfNorm = InfinityNorm( R );
            const double infResidual = RInfNorm / (AInfNorm*XInfNorm*epsilon*n);
            if( commRank == 0 )
            {
                if( details )
                {
                    Output("");
                    Output("||A||_oo       = ",AInfNorm);    
                    Output("||B||_oo       = ",BInfNorm);
                    Output("||X||_oo       = ",XInfNorm);
                    Output("||A X - B||_oo = ",RInfNorm);
                }
                Output
                ("||A X - B||_oo / (||A||_oo ||X||_oo eps n) = ",infResidual);
            }

            // Compute the relevant one norms and a relative residual
            const double AOneNorm = OneNorm( A );
            const double BOneNorm = OneNorm( B );
            const double XOneNorm = OneNorm( X );
            const double ROneNorm = OneNorm( R );
            const double oneResidual = ROneNorm / (AOneNorm*XOneNorm*epsilon*n);
            if( commRank == 0 )
            {
                if( details )
                {
                    Output("");
                    Output("||A||_1       = ",AOneNorm);    
                    Output("||B||_1       = ",BOneNorm);
                    Output("||X||_1       = ",XOneNorm);
                    Output("||A X - B||_1 = ",ROneNorm);
                }
                Output
                ("||A X - B||_1 / (||A||_1 ||X||_1 eps n) = ",oneResidual,"\n");
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
