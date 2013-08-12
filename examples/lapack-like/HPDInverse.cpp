/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/blas-like/level3/Hemm.hpp"
#include "elemental/lapack-like/Inverse.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
#include "elemental/matrices/Identity.hpp"
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
        const Int n = Input("--size","size of HPD matrix",100);
        const bool upper = Input("--upper","upper storage?",false);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        Grid g( mpi::COMM_WORLD );
        DistMatrix<C> A( g );
        HermitianUniformSpectrum( A, n, R(1), R(20) );

        if( print )
            Print( A, "A" );

        // Make a copy of A and then overwrite it with its inverse
        const UpperOrLower uplo = ( upper ? UPPER : LOWER );
        DistMatrix<C> invA( A );
        HPDInverse( uplo, invA );

        if( print )
        {
            MakeHermitian( uplo, invA );
            Print( invA, "inv(A)" );
        }

        // Form I - invA*A and print the relevant norms
        DistMatrix<C> E( g );
        Identity( E, n, n );
        Hemm( LEFT, uplo, C(-1), invA, A, C(1), E );

        const R frobNormA = HermitianFrobeniusNorm( uplo, A );
        const R frobNormInvA = HermitianFrobeniusNorm( uplo, invA );
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

