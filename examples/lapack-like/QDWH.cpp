/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental.hpp"
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

        Grid g( comm );
        DistMatrix<C> A( g ), Q( g ), P( g );
        Uniform( m, n, A );
        const R lowerBound = 1e-7;
        const R frobA = Norm( A, FROBENIUS_NORM );
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
        const int numItsQDWH = QDWH( Q, lowerBound, upperBound );
        Zeros( n, n, P );
        Gemm( ADJOINT, NORMAL, C(1), Q, A, C(0), P );

        // Check and report overall and orthogonality error
        DistMatrix<C> B( A );
        Gemm( NORMAL, NORMAL, C(-1), Q, P, C(1), B );
        const R frobQDWH = Norm( B, FROBENIUS_NORM );
        Identity( n, n, B );
        Herk( LOWER, NORMAL, C(1), Q, C(-1), B );
        const R frobQDWHOrthog = HermitianNorm( LOWER, B, FROBENIUS_NORM );
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
        const int numItsHalley = Halley( Q, upperBound );
        Zeros( n, n, P );
        Gemm( ADJOINT, NORMAL, C(1), Q, A, C(0), P );

        // Check and report the overall and orthogonality error
        B = A; 
        Gemm( NORMAL, NORMAL, C(-1), Q, P, C(1), B );
        const R frobHalley = Norm( B, FROBENIUS_NORM );
        Identity( n, n, B );
        Herk( LOWER, NORMAL, C(1), Q, C(-1), B );
        const R frobHalleyOrthog = HermitianNorm( LOWER, B, FROBENIUS_NORM );
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
