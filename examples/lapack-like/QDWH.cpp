/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/TwoUpperBound.hpp"
#include "elemental/lapack-like/Polar.hpp"
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

    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try 
    {
        const int m = Input("--height","height of matrix",100);
        const int n = Input("--width","width of matrix",100);
        ProcessInput();
        PrintInputReport();

        Grid g( comm );
        DistMatrix<C> A( g ), Q( g ), P( g );
        Uniform( m, n, A );
        const R lowerBound = 1e-7;
        const R frobA = FrobeniusNorm( A );
        const R upperBound = TwoNormUpperBound( A );
        if( g.Rank() == 0 )
        {
            std::cout << "ASSUMING 1 / ||inv(A)||_2 >= " << lowerBound << "\n"
                      << "||A||_F =  " << frobA << "\n"
                      << "||A||_2 <= " << upperBound << "\n" << std::endl;
        }

        // Compute the polar decomp of A using a QR-based Dynamically Weighted
        // Halley (QDWH) iteration
        Q = A;
        const int numItsQDWH = polar::QDWH( Q, lowerBound, upperBound );
        Zeros( n, n, P );
        Gemm( ADJOINT, NORMAL, C(1), Q, A, C(0), P );

        // Check and report overall and orthogonality error
        DistMatrix<C> B( A );
        Gemm( NORMAL, NORMAL, C(-1), Q, P, C(1), B );
        const R frobQDWH = FrobeniusNorm( B );
        Identity( n, n, B );
        Herk( LOWER, NORMAL, C(1), Q, C(-1), B );
        const R frobQDWHOrthog = HermitianFrobeniusNorm( LOWER, B );
        if( g.Rank() == 0 )
        {
            std::cout << numItsQDWH << " iterations of QDWH\n"
                      << "||A - QP||_F / ||A||_F = " 
                      << frobQDWH/frobA << "\n"
                      << "||I - QQ^H||_F / ||A||_F = " 
                      << frobQDWHOrthog/frobA << "\n"
                      << std::endl;
        }

        // Compute the polar decomp of A using a standard QR-based Halley
        // iteration
        Q = A;
        const int numItsHalley = polar::Halley( Q, upperBound );
        Zeros( n, n, P );
        Gemm( ADJOINT, NORMAL, C(1), Q, A, C(0), P );

        // Check and report the overall and orthogonality error
        B = A; 
        Gemm( NORMAL, NORMAL, C(-1), Q, P, C(1), B );
        const R frobHalley = FrobeniusNorm( B );
        Identity( n, n, B );
        Herk( LOWER, NORMAL, C(1), Q, C(-1), B );
        const R frobHalleyOrthog = HermitianFrobeniusNorm( LOWER, B );
        if( g.Rank() == 0 )
        {
            std::cout << numItsHalley << " iterations of Halley\n"
                      << "||A - QP||_F / ||A||_F = " 
                      << frobHalley/frobA << "\n"
                      << "||I - QQ^H||_F / ||A||_F = "
                      << frobHalleyOrthog/frobA << "\n"
                      << std::endl;
        }
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught exception with message: "
           << e.what() << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
