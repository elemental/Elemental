/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_IDENTITY_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const Int n = Input("--size","size of matrices to test",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( mpi::WorldRank() == 0 )
        {
            const Int commSize = mpi::Size( comm );
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
                Print( X, "Built-in identity" );
        }

        // Local buffers
        {
            // Allocate local data
            const Int gridHeight = grid.Height();
            const Int gridWidth = grid.Width();
            const Int gridRow = grid.Row();
            const Int gridCol = grid.Col();
            const Int localHeight = Length( n, gridRow, gridHeight );
            const Int localWidth = Length( n, gridCol, gridWidth );
            std::vector<double> localData( localHeight*localWidth );

            // Fill local data for identity
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                // Form global column index from local column index
                const Int j = gridCol + jLoc*gridWidth;
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                {
                    // Form global row index from local row index
                    const Int i = gridRow + iLoc*gridHeight;     
                    // If diagonal entry, set to one, otherwise zero
                    if( i == j )
                        localData[iLoc+jLoc*localHeight] = 1.;
                    else
                        localData[iLoc+jLoc*localHeight] = 0.;
                }
            }

            DistMatrix<double> X(grid);
            X.Attach( n, n, grid, 0, 0, localData.data(), localHeight );
            if( print )
                Print( X, "Identity constructed from local buffers" );

            // Build another set of local buffers and attach it to X.
            // This time, make it all two's.
            std::vector<double> localTwos( localHeight*localWidth, 2 ); 
            X.Attach( n, n, grid, 0, 0, localTwos.data(), localHeight );
            if( print )
                Print( X, "After viewing local buffers of all two's" );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
