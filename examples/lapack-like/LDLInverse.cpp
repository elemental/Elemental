/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/MakeSymmetric.hpp"
#include "elemental/blas-like/level3/Hemm.hpp"
#include "elemental/blas-like/level3/Symm.hpp"
#include "elemental/blas-like/level3/Trdtrmm.hpp"
#include "elemental/lapack-like/LDL.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/TriangularInverse.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Wigner.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int n = Input("--size","size of matrix to factor",100);
        const double realMean = Input("--realMean","real mean",0.);
        const double imagMean = Input("--imagMean","imag mean",0.);
        const double stddev = Input("--stddev","standard dev.",1.);
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

        // Make a copy of A and then overwrite it with its inverse
        // WARNING: There is no pivoting here!
        DistMatrix<C> invA( A );
        if( conjugate )
            LDLH( invA );
        else
            LDLT( invA );
        TriangularInverse( LOWER, UNIT, invA );
        Trdtrmm( orientation, LOWER, invA );

        // Form I - invA*A and print the relevant norms
        DistMatrix<C> E;
        Identity( E, n, n );
        if( conjugate )
            Hemm( LEFT, LOWER, C(-1), invA, A, C(1), E );
        else
            Symm( LEFT, LOWER, C(-1), invA, A, C(1), E );

        const Real frobNormA = FrobeniusNorm( A );
        const Real frobNormInvA = 
            ( conjugate ? HermitianFrobeniusNorm( LOWER, invA )
                        : SymmetricFrobeniusNorm( LOWER, invA ) );
        const Real frobNormError = FrobeniusNorm( E );
        if( mpi::WorldRank() == 0 )
        {
            std::cout << "|| A          ||_F = " << frobNormA << "\n"
                      << "|| invA       ||_F = " << frobNormInvA << "\n"
                      << "|| I - invA A ||_F = " << frobNormError << "\n"
                      << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
