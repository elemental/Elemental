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

void Usage()
{
    cout << "HPDInverse <n>\n"
         << "  <n>: size of random matrix to test HPDInverse on\n"
         << endl;
}

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    if( argc < 2 )
    {
        if( commRank == 0 )
            Usage();
        Finalize();
        return 0;
    }
    const int n = atoi( argv[1] );

    try 
    {
        Grid g( comm );
        DistMatrix<C> A( g );
        HermitianUniformSpectrum( n, A, R(1), R(20) );

        // Make a copy of A and then overwrite it with its inverse
        DistMatrix<C> invA( A );
        HPDInverse( LOWER, invA );

        // Form I - invA*A and print the relevant norms
        DistMatrix<C> E( g );
        Identity( n, n, E );
        Hemm( LEFT, LOWER, C(-1), invA, A, C(1), E );

        const R frobNormA = HermitianNorm( LOWER, A, FROBENIUS_NORM );
        const R frobNormInvA = HermitianNorm( LOWER, invA, FROBENIUS_NORM );
        const R frobNormError = Norm( E, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            std::cout << "|| A          ||_F = " << frobNormA << "\n"
                      << "|| invA       ||_F = " << frobNormInvA << "\n"
                      << "|| I - invA A ||_F = " << frobNormError << "\n"
                      << std::endl;
        }
    }
    catch( exception& e )
    {
        cerr << "Process " << commRank << " caught exception with message: "
             << e.what() << endl;
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}

