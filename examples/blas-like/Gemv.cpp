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
    cout << "Gemv <adjoint?> <m> <n> <print?>\n"
         << "  <adjoint?>: form A^H x if nonzero, A x otherwise\n"
         << "  <m>: height of matrix\n"
         << "  <n>: width of matrix\n"
         << "  <print?>: print matrices for verifying correctness\n"
         << endl;
}

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    if( argc < 5 )
    {
        if( commRank == 0 )
            Usage();
        Finalize();
        return 0;
    }
    const Orientation orientation = ( atoi(argv[1]) ? ADJOINT : NORMAL );
    const int m = atoi( argv[2] );
    const int n = atoi( argv[3] );
    const bool print = atoi( argv[4] );

    try 
    {
        Grid g( comm );
        DistMatrix<C> A( g );
        Uniform( m, n, A );

        // Draw the entries of the original x and y from uniform distributions 
        // over the complex unit ball
        DistMatrix<C> x( g ), y( g );
        if( orientation == NORMAL )
        {
            Uniform( n, 1, x );
            Uniform( m, 1, y );
        }
        else
        {
            Uniform( m, 1, x );
            Uniform( n, 1, y );
        }

        if( print )
        {
            A.Print("A");
            x.Print("x");
            y.Print("y");
        }

        // Run the matrix-vector product
        Gemv( orientation, (C)3, A, x, (C)4, y );

        if( print )
        {
            if( orientation == NORMAL )
                y.Print("y := 3 A x + 4 y");
            else
                y.Print("y := 3 A^H x + 4 y");
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

