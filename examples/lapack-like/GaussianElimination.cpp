/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace elem;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    try 
    {
        const int n = Input("--size","size of matrix",100);
        const int numRhs = Input("--numRhs","# of right-hand sides",1); 
        const int blocksize = Input("--blocksize","algorithmic blocksize",64);
        int gridHeight = Input("--gridHeight","grid height",0);
        const bool details = Input("--details","print norm details?",false);
        ProcessInput();
        PrintInputReport();

        // Set the algorithmic blocksize
        SetBlocksize( blocksize );

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( commSize );
        const int gridWidth = commSize / gridHeight;
        Grid grid( comm, gridHeight, gridWidth );

        // Set up random A and B, then make the copies X := B and ACopy := A
        DistMatrix<double> A(grid), B(grid), ACopy(grid), X(grid);
        for( int test=0; test<3; ++test )
        {
            Uniform( n, n,      A );
            Uniform( n, numRhs, B );
            ACopy = A;
            X = B;

            // Perform the LU factorization and simultaneous solve
            if( commRank == 0 )
            {
                std::cout << "Starting GaussianElimination...";
                std::cout.flush();
            }
            mpi::Barrier( comm );
            double startTime = mpi::Time();
            GaussianElimination( A, X );
            mpi::Barrier( comm );
            double stopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << stopTime-startTime << " seconds." << std::endl;

            // Form R := A X - B
            DistMatrix<double> R( B );
            Gemm( NORMAL, NORMAL, (double)1, ACopy, X, (double)-1, R );

            // Compute the relevant Frobenius norms and a relative residual
            const double epsilon = lapack::MachineEpsilon<double>();
            const double AFrobNorm = Norm( ACopy, FROBENIUS_NORM );
            const double BFrobNorm = Norm( B,     FROBENIUS_NORM );
            const double XFrobNorm = Norm( X,     FROBENIUS_NORM );
            const double RFrobNorm = Norm( R,     FROBENIUS_NORM );
            const double frobResidual = 
                RFrobNorm / (AFrobNorm*XFrobNorm*epsilon*n);
            if( commRank == 0 )
            {
                if( details )
                    std::cout << "||A||_F       = " << AFrobNorm << "\n"
                              << "||B||_F       = " << BFrobNorm << "\n"
                              << "||X||_F       = " << XFrobNorm << "\n"
                              << "||A X - B||_F = " << RFrobNorm << "\n";
                std::cout << "||A X - B||_F / (||A||_F ||X||_F epsilon n) = " 
                          << frobResidual << "\n";
            }

            // Compute the relevant infinity norms and a relative residual
            const double AInfNorm = Norm( ACopy, INFINITY_NORM );
            const double BInfNorm = Norm( B,     INFINITY_NORM );
            const double XInfNorm = Norm( X,     INFINITY_NORM );
            const double RInfNorm = Norm( R,     INFINITY_NORM );
            const double infResidual = RInfNorm / (AInfNorm*XInfNorm*epsilon*n);
            if( commRank == 0 )
            {
                if( details )
                    std::cout << "\n"
                              << "||A||_oo       = " << AInfNorm << "\n"
                              << "||B||_oo       = " << BInfNorm << "\n"
                              << "||X||_oo       = " << XInfNorm << "\n"
                              << "||A X - B||_oo = " << RInfNorm << "\n";
                std::cout << "||A X - B||_oo / (||A||_oo ||X||_oo epsilon n) = "
                          << infResidual << "\n";
            }

            // Compute the relevant one norms and a relative residual
            const double AOneNorm = Norm( ACopy, ONE_NORM );
            const double BOneNorm = Norm( B,     ONE_NORM );
            const double XOneNorm = Norm( X,     ONE_NORM );
            const double ROneNorm = Norm( R,     ONE_NORM );
            const double oneResidual = ROneNorm / (AOneNorm*XOneNorm*epsilon*n);
            if( commRank == 0 )
            {
                if( details )
                    std::cout << "\n"
                              << "||A||_1       = " << AOneNorm << "\n"
                              << "||B||_1       = " << BOneNorm << "\n"
                              << "||X||_1       = " << XOneNorm << "\n"
                              << "||A X - B||_1 = " << ROneNorm << "\n";
                std::cout << "||A X - B||_1 / (||A||_1 ||X||_1 epsilon n) = " 
                          << oneResidual << "\n" << std::endl;
            }
        }
    }
    catch( ArgException& e )
    {
        // There is no reason to do anything
    }
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught exception: " << e.what()
           << std::endl;
        std::cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
