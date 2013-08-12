/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Hemm.hpp"
#include "elemental/blas-like/level3/Symm.hpp"
#include "elemental/blas-like/level3/Trdtrmm.hpp"
#include "elemental/lapack-like/LDL.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/TriangularInverse.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int n = Input("--size","size of matrix to factor",100);
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        ProcessInput();
        PrintInputReport();

        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        Grid g( mpi::COMM_WORLD );
        DistMatrix<C> A( g );

        if( conjugate )
        {
            HermitianUniformSpectrum( A, n, -30, -20 );
        }
        else
        {
            Uniform( A, n, n );
            DistMatrix<C> ATrans( g );
            Transpose( A, ATrans );
            Axpy( C(1), ATrans, A );
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
        DistMatrix<C> E( g );
        Identity( E, n, n );
        if( conjugate )
            Hemm( LEFT, LOWER, C(-1), invA, A, C(1), E );
        else
            Symm( LEFT, LOWER, C(-1), invA, A, C(1), E );

        const R frobNormA = FrobeniusNorm( A );
        const R frobNormInvA = 
            ( conjugate ? HermitianFrobeniusNorm( LOWER, invA )
                        : SymmetricFrobeniusNorm( LOWER, invA ) );
        const R frobNormError = FrobeniusNorm( E );
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
