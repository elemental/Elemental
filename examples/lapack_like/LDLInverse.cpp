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
        const Int n = Input("--size","size of matrix to factor",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Real realMean = Input("--realMean","real mean",Real(0));
        const Real imagMean = Input("--imagMean","imag mean",Real(0));
        const Real stddev = Input("--stddev","standard dev.",Real(1));
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );

        C mean( realMean, imagMean );
        DistMatrix<C> A;
        if( conjugate )
        {
            Wigner( A, n, mean, stddev );
            //HermitianUniformSpectrum( A, n, 1, 2 );
        }
        else
        {
            Gaussian( A, n, n, mean, stddev );
            MakeSymmetric( LOWER, A );
        }

        Timer timer;
        // Make a copy of A and then overwrite it with its inverse
        DistMatrix<C> invA( A );
        if( mpi::Rank() == 0 )
            timer.Start();
        SymmetricInverse( LOWER, invA, conjugate );
        if( mpi::Rank() == 0 )
            timer.Stop();

        // Form I - invA*A and print the relevant norms
        DistMatrix<C> E;
        Identity( E, n, n );
        Symm( LEFT, LOWER, C(-1), invA, A, C(1), E, conjugate );

        const Real frobNormA = FrobeniusNorm( A );
        const Real frobNormInvA = SymmetricFrobeniusNorm( LOWER, invA );
        const Real frobNormError = FrobeniusNorm( E );
        if( mpi::Rank() == 0 )
        {
            Output("LDLInverse time: ",timer.Total()," secs");
            Output
            ("|| A          ||_F = ",frobNormA,"\n",
             "|| inv(A)     ||_F = ",frobNormInvA,"\n",
             "|| I - inv(A) ||_F = ",frobNormError,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
