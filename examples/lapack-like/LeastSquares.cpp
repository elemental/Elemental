/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/LeastSquares.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/matrices/Uniform.hpp"
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
        const bool adjoint = Input("--adjoint","adjoint solve?",false);
        const int m = Input("--height","height of matrix",100);
        const int n = Input("--width","width of matrix",100);
        const int numRhs = Input("--numRhs","# of right-hand sides",1);
        const int blocksize = Input("--blocksize","algorithmic blocksize",64);
        int gridHeight = Input("--gridHeight","grid height",0);
        ProcessInput();
        PrintInputReport();

        const Orientation orientation = ( adjoint ? ADJOINT : NORMAL );

        // Set the algorithmic blocksize
        SetBlocksize( blocksize );

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( commSize );
        const int gridWidth = commSize / gridHeight;
        Grid grid( comm, gridHeight, gridWidth );

        // Set up random A and B, then make the copies X := B and ACopy := A
        typedef Complex<double> F;
        DistMatrix<F> A(grid), B(grid), ACopy(grid), X(grid), Z(grid);
        for( int test=0; test<3; ++test )
        {
            const int k = ( orientation==NORMAL ? m : n );
            const int N = ( orientation==NORMAL ? n : m );
            Uniform( m, n, A );
            Zeros( k, numRhs, B );
            ACopy = A;

            // Form B in the range of op(A)
            Uniform( N, numRhs, Z );
            Gemm( orientation, NORMAL, F(1), A, Z, F(0), B );

            // Perform the QR/LQ factorization and solve
            if( commRank == 0 )
            {
                std::cout << "Starting LeastSquares...";
                std::cout.flush();
            }
            mpi::Barrier( comm );
            double startTime = mpi::Time();
            LeastSquares( orientation, A, B, X );
            mpi::Barrier( comm );
            double stopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << stopTime-startTime << " seconds." << std::endl;

            // Form R := op(A) X - B
            DistMatrix<F> R( B );
            Gemm( orientation, NORMAL, F(1), ACopy, X, F(-1), R );

            // Compute the relevant Frobenius norms and a relative residual
            const double epsilon = lapack::MachineEpsilon<double>();
            const double AFrobNorm = FrobeniusNorm( ACopy );
            const double BFrobNorm = FrobeniusNorm( B );
            const double XFrobNorm = FrobeniusNorm( X );
            const double RFrobNorm = FrobeniusNorm( R );
            const double frobResidual = 
                RFrobNorm / (AFrobNorm*XFrobNorm*epsilon*n);
            if( commRank == 0 )
                std::cout << "||A||_F       = " << AFrobNorm << "\n"
                          << "||B||_F       = " << BFrobNorm << "\n"
                          << "||X||_F       = " << XFrobNorm << "\n"
                          << "||A X - B||_F = " << RFrobNorm << "\n"
                          << "||A X - B||_F / (||A||_F ||X||_F epsilon n) = " 
                          << frobResidual << "\n" << std::endl;

            // Compute the relevant infinity norms and a relative residual
            const double AInfNorm = InfinityNorm( ACopy );
            const double BInfNorm = InfinityNorm( B );
            const double XInfNorm = InfinityNorm( X );
            const double RInfNorm = InfinityNorm( R );
            const double infResidual = RInfNorm / (AInfNorm*XInfNorm*epsilon*n);
            if( commRank == 0 )
                std::cout << "||A||_oo       = " << AInfNorm << "\n"
                          << "||B||_oo       = " << BInfNorm << "\n"
                          << "||X||_oo       = " << XInfNorm << "\n"
                          << "||A X - B||_oo = " << RInfNorm << "\n"
                          << "||A X - B||_oo / (||A||_oo ||X||_oo epsilon n) = "
                          << infResidual << "\n" << std::endl;

            // Compute the relevant one norms and a relative residual
            const double AOneNorm = OneNorm( ACopy );
            const double BOneNorm = OneNorm( B );
            const double XOneNorm = OneNorm( X );
            const double ROneNorm = OneNorm( R );
            const double oneResidual = ROneNorm / (AOneNorm*XOneNorm*epsilon*n);
            if( commRank == 0 )
                std::cout << "||A||_1       = " << AOneNorm << "\n"
                          << "||B||_1       = " << BOneNorm << "\n"
                          << "||X||_1       = " << XOneNorm << "\n"
                          << "||A X - B||_1 = " << ROneNorm << "\n"
                          << "||A X - B||_1 / (||A||_1 ||X||_1 epsilon n) = " 
                          << oneResidual << "\n" << std::endl;
            
            if( commRank == 0 )
                std::cout << std::endl;
        }
    }
    catch( ArgException& e )
    {
        // There is nothing to do
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
