/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
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
        const Int n = Input("--size","size of matrix to factor",100);
        const Real realMean = Input("--realMean","real mean",Real(0)); 
        const Real imagMean = Input("--imagMean","imag mean",Real(0));
        const Real stddev = Input("--stddev","standard dev.",Real(1));
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        ProcessInput();
        PrintInputReport();

        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        C mean( realMean, imagMean );
        DistMatrix<C> A;
        if( conjugate )
        {
            Wigner( A, n, mean, stddev );
        }
        else
        {
            Gaussian( A, n, n, mean, stddev );
            MakeSymmetric( LOWER, A );
        }

        // Make a copy of A and then overwrite it with its LDL factorization
        // WARNING: There is no pivoting here!
        Timer timer;
        DistMatrix<C> factA( A );
        if( mpi::Rank() == 0 )
            timer.Start();
        LDL( factA, conjugate );
        if( mpi::Rank() == 0 )
            timer.Stop();
        auto d = GetDiagonal(factA);

        DistMatrix<C> L( factA );
        MakeTrapezoidal( LOWER, L );
        FillDiagonal( L, C(1) );

        DistMatrix<C> LD( L );
        DiagonalScale( RIGHT, NORMAL, d, LD );
        Gemm( NORMAL, orientation, C(-1), LD, L, C(1), A );
        const Real frobNormError = FrobeniusNorm( A );
        if( mpi::Rank() == 0 )
        {
            Output("LDL time: ",timer.Total()," secs");
            Output("|| A - L D L^[T/H] ||_F = ",frobNormError,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
