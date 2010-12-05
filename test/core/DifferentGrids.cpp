/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include <cstdlib>
#include <ctime>
#include "elemental.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Run some tests for creating matrices with different grids.\n\n"
         << "  DifferentGrids <m> <n>\n\n"
         << "  m: height of matrices\n"
         << "  n: width of matrices\n" << endl;
}

int main( int argc, char* argv[] )
{
    int rank;
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 3 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }
    const int m = atoi(argv[1]);
    const int n = atoi(argv[2]);
#ifndef RELEASE
    if( rank == 0 )
    {
        cout << "==========================================\n"
             << " In debug mode! Performance will be poor! \n"
             << "==========================================" << endl;
    }
#endif
    try
    {
        int p;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Comm_size( comm, &p );

        std::vector<int> evenRanks((p+1)/2);
        std::vector<int> oddRanks(p/2);
        for( int i=0; i<evenRanks.size(); ++i )
            evenRanks[i] = 2*i;
        for( int i=0; i<oddRanks.size(); ++i )
            oddRanks[i] = 2*i + 1;

        MPI_Group group, evenGroup, oddGroup;
        MPI_Comm_group( comm, &group );
        MPI_Group_incl( group, evenRanks.size(), &evenRanks[0], &evenGroup );
        MPI_Group_incl( group, oddRanks.size(), &oddRanks[0], &oddGroup );

        if( rank == 0 )
        {
            std::cout << "Creating even grid...";
            std::cout.flush();
        }
        const Grid evenGrid( comm, evenGroup );
        if( rank == 0 )
            std::cout << "done." << std::endl;
        if( rank == 0 )
        {
            std::cout << "Creating odd grid...";
            std::cout.flush();
        }
        const Grid oddGrid( comm, oddGroup );
        if( rank == 0 )
            std::cout << "done." << std::endl;

        /*
        DistMatrix<double,MC,MR> AEven( m, n, evenGrid );
        DistMatrix<double,MC,MR> AOdd( m, n, oddGrid );

        AEven.Print("AEven");
        AOdd.Print("AOdd");
        */
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

