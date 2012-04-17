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
#include <cstdlib>
#include <ctime>
#include "elemental.hpp"
using namespace elem;

void Usage()
{
    std::cout 
        << "Run some tests for creating matrices with different grids.\n\n"
        << "  DifferentGrids <m> <n>\n\n"
        << "  m: height of matrices\n"
        << "  n: width of matrices\n" << std::endl;
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    
    if( argc < 3 )
    {
        if( commRank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try
    {
        int argNum = 0;
        const int m = atoi(argv[++argNum]);
        const int n = atoi(argv[++argNum]);
#ifndef RELEASE
        if( commRank == 0 )
        {
            std::cout 
                << "==========================================\n"
                << " In debug mode! Performance will be poor! \n"
                << "==========================================" << std::endl;
        }
#endif

        // Drop down to a square grid, change the matrix, and redistribute back
        const int commSqrt = 
            static_cast<int>(sqrt(static_cast<double>(commSize)));

        std::vector<int> sqrtRanks(commSqrt*commSqrt);
        for( int i=0; i<commSqrt*commSqrt; ++i )
            sqrtRanks[i] = i;

        mpi::Group group, sqrtGroup;
        
        mpi::CommGroup( comm, group );
        mpi::GroupIncl( group, sqrtRanks.size(), &sqrtRanks[0], sqrtGroup );

        const Grid grid( comm );
        const Grid sqrtGrid( comm, sqrtGroup );

        DistMatrix<double> A(grid), ASqrt(sqrtGrid);

        Identity( m, n, A );
        A.Print("A");

        ASqrt = A;
        ASqrt.Print("ASqrt := A");

        Scal( 2.0, ASqrt );
        ASqrt.Print("ASqrt := 2 ASqrt");

        A = ASqrt;
        A.Print("A := ASqrt");
    }
    catch( std::exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        std::cerr << "Process " << commRank << " caught error message:\n"
                  << e.what() << std::endl;
    }
    Finalize();
    return 0;
}

