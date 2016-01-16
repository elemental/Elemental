/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

typedef double Real;
typedef Complex<Real> F;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
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
        const bool elemental = Input("--elemental","test Elemental?",true);
        const bool error = Input("--error","test Elemental error?",true);
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
        DistMatrix<F> A(grid), B(grid), X(grid);

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
                    Output("Starting ScaLAPACK linear solve");
                DistMatrix<F,MC,MR,BLOCK> ABlock( A ), BBlock( B );
                mpi::Barrier( comm );
                if( commRank == 0 )
                    timer.Start();
                LinearSolve( ABlock, BBlock );
                if( commRank == 0 )
                    Output(timer.Stop()," seconds");
                if( error && !elemental )
                    X = BBlock;
            }

            // Perform the LU factorization and simultaneous solve
            if( elemental )
            {
                if( commRank == 0 )
                    Output("Starting Elemental linear solve");
                mpi::Barrier( comm );
                if( commRank == 0 )
                    timer.Start();
                LinearSolve( A, X );
                mpi::Barrier( comm );
                if( commRank == 0 )
                    Output(timer.Stop()," seconds");
            }

            if( error && (elemental || scalapack) )
            { 
                // Form R := A X - B
                auto R( B );
                Gemm( NORMAL, NORMAL, F(1), A, X, F(-1), R );

                // Compute infinity norms and a relative residual
                const Real eps = limits::Epsilon<Real>();
                const Real AInfNorm = InfinityNorm( A );
                const Real BInfNorm = InfinityNorm( B );
                const Real XInfNorm = InfinityNorm( X );
                const Real RInfNorm = InfinityNorm( R );
                const Real infResidual = RInfNorm / (AInfNorm*XInfNorm*eps*n);
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
                    ("||A X - B||_oo / (||A||_oo ||X||_oo eps n) = ",
                     infResidual);
                }

                // Compute one norms and a relative residual
                const Real AOneNorm = OneNorm( A );
                const Real BOneNorm = OneNorm( B );
                const Real XOneNorm = OneNorm( X );
                const Real ROneNorm = OneNorm( R );
                const Real oneResidual = ROneNorm / (AOneNorm*XOneNorm*eps*n);
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
                    ("||A X - B||_1 / (||A||_1 ||X||_1 eps n) = ",
                     oneResidual,"\n");
                }
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
