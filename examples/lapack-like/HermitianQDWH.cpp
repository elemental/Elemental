/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Polar.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/TwoUpperBound.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
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
        const Int n = Input("--size","size of Hermitian matrix",100);
        ProcessInput();
        PrintInputReport();

        Grid g( mpi::COMM_WORLD );
        DistMatrix<C> A( g ), Q( g ), P( g );
        HermitianUniformSpectrum( A, n, 1, 2 );
        const R lowerBound = 1.0;
        const R frobA = FrobeniusNorm( A );
        const R upperBound = TwoNormUpperBound( A );
        if( mpi::WorldRank() == 0 )
        {
            std::cout << "ASSUMING 1 / ||inv(A)||_2 >= " << lowerBound << "\n"
                      << "||A||_F =  " << frobA << "\n"
                      << "||A||_2 <= " << upperBound << "\n" << std::endl;
        }

        // Compute the polar decomp of A using a QR-based Dynamically Weighted
        // Halley (QDWH) iteration
        Q = A;
        const Int numItsQDWH = 
            hermitian_polar::QDWH( LOWER, Q, lowerBound, upperBound );
        Zeros( P, n, n );
        Gemm( ADJOINT, NORMAL, C(1), Q, A, C(0), P );

        // Check and report overall and orthogonality error
        DistMatrix<C> B( A );
        Gemm( NORMAL, NORMAL, C(-1), Q, P, C(1), B );
        const R frobQDWH = FrobeniusNorm( B );
        Identity( B, n, n );
        Herk( LOWER, NORMAL, C(1), Q, C(-1), B );
        const R frobQDWHOrthog = HermitianFrobeniusNorm( LOWER, B );
        if( mpi::WorldRank() == 0 )
        {
            std::cout << numItsQDWH << " iterations of QDWH\n"
                      << "||A - QP||_F / ||A||_F = " 
                      << frobQDWH/frobA << "\n"
                      << "||I - QQ^H||_F / ||A||_F = " 
                      << frobQDWHOrthog/frobA << "\n"
                      << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
