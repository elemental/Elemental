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
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;

void Usage()
{
    cout << "AxpyInterface <n>\n"
         << "\n"
         << "  <n>: size of matrix to test.\n"
         << std::endl;
}

int
main( int argc, char* argv[] )
{
    Init( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    int rank = mpi::CommRank( comm );

    if( argc < 2 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try 
    {
        const int n = atoi(argv[1]);

        Grid g( comm );

        DistMatrix<double,MC,MR> A( n, n, g );
        A.SetToZero();

        AxpyInterface<double> interface;
        interface.Attach( A );
        Matrix<double> X( n, n );
        const int rank = A.Grid().VCRank();
        for( int j=0; j<n; ++j )
            for( int i=0; i<n; ++i )
                X.Set(i,j,rank);
        interface.Axpy( 1.0, X, 0, 0 );
        interface.Detach();

        // Ensure that our local matrix is the sum of all the ranks
        const int p = mpi::CommSize( mpi::COMM_WORLD );
        const double sumOfRanks = ((p-1)*(p-1)+(p-1))/2;
        // Check that our local matrix is equal to sumOfRanks everywhere
        double myMaxError = 0.;
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
            for( int iLocal=0; iLocal<A.LocalHeight(); ++iLocal )
                myMaxError = 
                    std::max( myMaxError, 
                              Abs(sumOfRanks-A.GetLocalEntry(iLocal,jLocal)) );
        double maxError; 
        mpi::AllReduce
        ( &myMaxError, &maxError, 1, mpi::SUM, g.VCComm() );

        if( rank == 0 )
            std::cout << "max error = " << maxError << std::endl;
        if( maxError > 0.000001 )
            A.Print("A");
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

