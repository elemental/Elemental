/*
   Copyright (c) 2009-2011, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
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
using namespace elemental;

void Usage()
{
    cout << "AxpyInterface\n" << std::endl;
}

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );
    const int p = mpi::CommSize( comm );

    try 
    {
        const int m = 3*p;
        const int n = 2*p;

        Grid g( comm );

        for( int k=0; k<50; ++k )
        {
            if( rank == 0 )
                std::cout << "Iteration " << k << std::endl;

            DistMatrix<double,MC,MR> A( m, n, g );
            A.SetToZero();

            AxpyInterface<double> interface;
            interface.Attach( LOCAL_TO_GLOBAL, A );
            Matrix<double> X( p, 1 );
            for( int j=0; j<X.Width(); ++j )
                for( int i=0; i<p; ++i )
                    X.Set(i,j,rank+1);
            for( int i=0; i<5; ++i )
            {
                interface.Axpy( 2, X, 2*rank, rank );
                interface.Axpy( 2, X, 2*rank, rank+1 );
            }
            interface.Detach();

            A.Print("A");

            interface.Attach( GLOBAL_TO_LOCAL, A );
            Matrix<double> Y;
            if( rank == 0 )
            {
                Y.ResizeTo( m, n );
                Y.SetToZero();
                interface.Axpy( 1.0, Y, 0, 0 );
            }
            interface.Detach();

            if( rank == 0 )
                Y.Print( "Copy of global matrix on root process:" );

            // TODO: Check to ensure that the result is correct
        }
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:\n"
             << e.what() << endl;
    }

    Finalize();
    return 0;
}

