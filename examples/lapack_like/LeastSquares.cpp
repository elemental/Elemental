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

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try 
    {
        const char trans = Input("--trans","orientation",'N');
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int numRhs = Input("--numRhs","# of right-hand sides",1);
        const Int blocksize = Input("--blocksize","algorithmic blocksize",64);
        Int gridHeight = Input("--gridHeight","grid height",0);
        ProcessInput();
        PrintInputReport();

        const Orientation orientation = CharToOrientation( trans );

        // Set the algorithmic blocksize
        SetBlocksize( blocksize );

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( commSize );
        Grid grid( comm, gridHeight );

        // Set up random A and B, then make the copy X := B
        DistMatrix<F> A(grid), B(grid), X(grid), Z(grid);
        for( Int test=0; test<3; ++test )
        {
            const Int k = ( orientation==NORMAL ? m : n );
            const Int N = ( orientation==NORMAL ? n : m );
            Uniform( A, m, n );
            Zeros( B, k, numRhs );

            // Form B in the range of op(A)
            Uniform( Z, N, numRhs );
            Gemm( orientation, NORMAL, F(1), A, Z, F(0), B );

            // Perform the QR/LQ factorization and solve
            if( commRank == 0 )
                Output("Starting LeastSquares...");
            mpi::Barrier();
            double startTime = mpi::Time();
            LeastSquares( orientation, A, B, X );
            mpi::Barrier();
            double stopTime = mpi::Time();
            if( commRank == 0 )
                Output(stopTime-startTime," seconds.");

            // Form R := op(A) X - B
            DistMatrix<F> R( B );
            Gemm( orientation, NORMAL, F(1), A, X, F(-1), R );

            // Compute the relevant Frobenius norms and a relative residual
            const Real eps = limits::Epsilon<Real>();
            const Real AFrobNorm = FrobeniusNorm( A );
            const Real BFrobNorm = FrobeniusNorm( B );
            const Real XFrobNorm = FrobeniusNorm( X );
            const Real RFrobNorm = FrobeniusNorm( R );
            const Real frobResidual = RFrobNorm / (AFrobNorm*XFrobNorm*eps*n);
            if( commRank == 0 )
                Output
                ("|| A ||_F       = ",AFrobNorm,"\n",
                 "|| B ||_F       = ",BFrobNorm,"\n",
                 "|| X ||_F       = ",XFrobNorm,"\n",
                 "|| A X - B ||_F = ",RFrobNorm,"\n",
                 "|| op(A) X - B ||_F / (|| A ||_F || X ||_F eps n) = ",
                 frobResidual,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
