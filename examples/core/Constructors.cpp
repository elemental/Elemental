/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/matrices/Identity.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const int n = Input("--size","size of matrices to test",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( mpi::WorldRank() == 0 )
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
            Identity( X, n, n );
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
            const int localHeight = Length( n, gridRow, gridHeight );
            const int localWidth = Length( n, gridCol, gridWidth );
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
            X.Attach( n, n, 0, 0, &localTwos[0], localHeight, grid );
            if( print )
                X.Print("After viewing local buffers of all two's");
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
