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
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n = Input("--size","size of matrices to test",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();

        if( commRank == 0 )
        {
            const int commSize = mpi::CommSize( comm );
            std::cout << "Will create matrices distributed over " 
                      << commSize << " process(es) in various ways" 
                      << std::endl;
        }

        // Built-in
        const Grid grid( comm );
        {
            DistMatrix<double> X(grid);
            Identity( n, n, X );
            if( print )
                X.Print("Built-in identity");
        }

        // Local buffers
        {
            // Allocate local data
            const int gridHeight = grid.Height();
            const int gridWidth = grid.Width();
            const int gridRow = grid.Row();
            const int gridCol = grid.Col();
            const int localHeight = LocalLength( n, gridRow, gridHeight );
            const int localWidth = LocalLength( n, gridCol, gridWidth );
            std::vector<double> localData( localHeight*localWidth );

            // Fill local data for identity
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                // Form global column index from local column index
                const int j = gridCol + jLocal*gridWidth;
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                {
                    // Form global row index from local row index
                    const int i = gridRow + iLocal*gridHeight;     

                    // If diagonal entry, set to one, otherwise zero
                    if( i == j )
                        localData[iLocal+jLocal*localHeight] = 1.;
                    else
                        localData[iLocal+jLocal*localHeight] = 0.;
                }
            }

            DistMatrix<double> 
                X( n, n, 0, 0, &localData[0], localHeight, grid );
            if( print )
                X.Print("Identity constructed from local buffers");

            // Build another set of local buffers and attach it to X.
            // This time, make it all two's.
            std::vector<double> localTwos( localHeight*localWidth, 2 ); 
            X.View( n, n, 0, 0, &localTwos[0], localHeight, grid );
            if( print )
                X.Print("After viewing local buffers of all two's");
        }
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught error message:\n"
           << e.what() << std::endl;
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
