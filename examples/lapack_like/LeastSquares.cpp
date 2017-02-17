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
    const int commRank = El::mpi::Rank( comm );
    const int commSize = El::mpi::Size( comm );

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const char trans = El::Input("--trans","orientation",'N');
        const El::Int m = El::Input("--height","height of matrix",100);
        const El::Int n = El::Input("--width","width of matrix",100);
        const El::Int numRhs = El::Input("--numRhs","# of right-hand sides",1);
        const El::Int blocksize =
          El::Input("--blocksize","algorithmic blocksize",64);
        El::Int gridHeight = El::Input("--gridHeight","grid height",0);
        El::ProcessInput();
        El::PrintInputReport();

        const El::Orientation orientation = El::CharToOrientation( trans );

        // Set the algorithmic blocksize
        El::SetBlocksize( blocksize );

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = El::Grid::DefaultHeight( commSize );
        const El::Grid grid( comm, gridHeight );

        // Set up random A and B, then make the copy X := B
        El::DistMatrix<Scalar> A(grid), B(grid), X(grid), Z(grid);
        for( El::Int test=0; test<3; ++test )
        {
            const El::Int k = orientation==El::NORMAL ? m : n;
            const El::Int N = orientation==El::NORMAL ? n : m;
            El::Uniform( A, m, n );
            El::Zeros( B, k, numRhs );

            // Form B in the range of op(A)
            El::Uniform( Z, N, numRhs );
            El::Gemm( orientation, El::NORMAL, Scalar(1), A, Z, Scalar(0), B );

            // Perform the QR/LQ factorization and solve
            if( commRank == 0 )
                El::Output("Starting LeastSquares...");
            El::mpi::Barrier( comm );
            double startTime = El::mpi::Time();
            El::LeastSquares( orientation, A, B, X );
            El::mpi::Barrier();
            double stopTime = El::mpi::Time();
            if( commRank == 0 )
                El::Output(stopTime-startTime," seconds.");

            // Form R := op(A) X - B
            El::DistMatrix<Scalar> R( B );
            El::Gemm( orientation, El::NORMAL, Scalar(1), A, X, Scalar(-1), R );

            // Compute the relevant Frobenius norms and a relative residual
            const Real eps = El::limits::Epsilon<Real>();
            const Real AFrobNorm = El::FrobeniusNorm( A );
            const Real BFrobNorm = El::FrobeniusNorm( B );
            const Real XFrobNorm = El::FrobeniusNorm( X );
            const Real RFrobNorm = El::FrobeniusNorm( R );
            const Real frobResidual = RFrobNorm / (AFrobNorm*XFrobNorm*eps*n);
            if( commRank == 0 )
                El::Output
                ("|| A ||_F       = ",AFrobNorm,"\n",
                 "|| B ||_F       = ",BFrobNorm,"\n",
                 "|| X ||_F       = ",XFrobNorm,"\n",
                 "|| A X - B ||_F = ",RFrobNorm,"\n",
                 "|| op(A) X - B ||_F / (|| A ||_F || X ||_F eps n) = ",
                 frobResidual,"\n");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
