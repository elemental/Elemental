/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Field>
void DemonstrateConstructors( El::Int n, const El::Grid& grid, bool print )
{
    if( grid.Rank() == 0 )
        El::Output
        ("Testing constructors of DistMatrix<",El::TypeName<Field>(),">");

    // Create an identity matrix after using the default constructor
    {
        El::DistMatrix<Field> X(grid);
        El::Identity( X, n, n );
        if( print )
            El::Print( X, "Used built-in identity" );
    }

    // Create an identity matrix by attaching pre-filled local buffers
    {
        // Allocate local data
        const El::Int gridHeight = grid.Height();
        const El::Int gridWidth = grid.Width();
        const El::Int gridRow = grid.Row();
        const El::Int gridCol = grid.Col();
        const El::Int localHeight = El::Length( n, gridRow, gridHeight );
        const El::Int localWidth = El::Length( n, gridCol, gridWidth );
        std::vector<Field> localData( localHeight*localWidth );

        // Fill local data for identity
        for( El::Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            // Form global column index from local column index
            const El::Int j = gridCol + jLoc*gridWidth;
            for( El::Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                // Form global row index from local row index
                const El::Int i = gridRow + iLoc*gridHeight;
                // If diagonal entry, set to one, otherwise zero
                if( i == j )
                    localData[iLoc+jLoc*localHeight] = 1;
                else
                    localData[iLoc+jLoc*localHeight] = 0;
            }
        }

        El::DistMatrix<Field> X(grid);
        const int colAlign = 0;
        const int rowAlign = 0;
        X.Attach
        ( n, n, grid, colAlign, rowAlign, localData.data(), localHeight );
        if( print )
            El::Print( X, "Identity constructed from local buffers" );
    }
}

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int n = El::Input("--size","size of matrices to test",100);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        // We will distribute our work using all processes
        El::mpi::Comm comm = El::mpi::COMM_WORLD;
        if( El::mpi::Rank(comm) == 0 )
            El::Output
            ("Creating matrices distributed over ",El::mpi::Size(comm),
             " process(es) in various ways");

        // Construct the default grid over COMM_WORLD
        const El::Grid grid( comm );

        DemonstrateConstructors<float>
        ( n, grid, print );
        DemonstrateConstructors<El::Complex<float>>
        ( n, grid, print );
        DemonstrateConstructors<double>
        ( n, grid, print );
        DemonstrateConstructors<El::Complex<double>>
        ( n, grid, print );
#ifdef EL_HAVE_QD
        DemonstrateConstructors<El::DoubleDouble>
        ( n, grid, print );
        DemonstrateConstructors<El::Complex<El::DoubleDouble>>
        ( n, grid, print );

        DemonstrateConstructors<El::QuadDouble>
        ( n, grid, print );
        DemonstrateConstructors<El::Complex<El::QuadDouble>>
        ( n, grid, print );
#endif
#ifdef EL_HAVE_QUAD
        DemonstrateConstructors<El::Quad>
        ( n, grid, print );
        DemonstrateConstructors<El::Complex<El::Quad>>
        ( n, grid, print );
#endif
#ifdef EL_HAVE_MPC
        DemonstrateConstructors<El::BigFloat>
        ( n, grid, print );
        DemonstrateConstructors<El::Complex<El::BigFloat>>
        ( n, grid, print );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
