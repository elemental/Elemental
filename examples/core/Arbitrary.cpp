/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Real>
void CorruptFourierSubmatrix
( El::Int n,
  El::Int numCorruptRows,
  El::Int numCorruptCols,
  const El::Grid& grid,
  bool print,
  bool display )
{
    typedef El::Complex<Real> Field;
    if( grid.Rank() == 0 )
        std::cout << "Testing with " << El::TypeName<Field>() << std::endl;

    El::DistMatrix<Field> A(grid);
    El::Fourier( A, n );
    if( display )
        El::Display( A, "Fourier Matrix" );
    if( print )
        El::Print( A, "Fourier matrix:" );

    // Get a consistent set of row and column indices (duplication is okay)
    std::vector<El::Int> rowInds(numCorruptRows), colInds(numCorruptCols);
    const int randomRoot = 0;
    if( grid.Rank() == randomRoot )
    {
        for( El::Int j=0; j<numCorruptRows; ++j )
            rowInds[j] = El::SampleUniform(El::Int(0),n);
        for( El::Int j=0; j<numCorruptCols; ++j )
            colInds[j] = El::SampleUniform(El::Int(0),n);
    }
    El::mpi::Broadcast
    ( rowInds.data(), numCorruptRows, randomRoot, grid.Comm() );
    El::mpi::Broadcast
    ( colInds.data(), numCorruptCols, randomRoot, grid.Comm() );
    if( grid.Rank() == randomRoot && print )
    {
        std::cout << "rowInds: \n";
        for( El::Int j=0; j<numCorruptRows; ++j )
            std::cout << rowInds[j] << "\n";
        std::cout << "\n";
        std::cout << "colInds: \n";
        for( El::Int j=0; j<numCorruptCols; ++j )
            std::cout << colInds[j] << "\n";
        std::cout << std::endl;
    }

    auto ASub =  A( rowInds, colInds );
    if( display )
        El::Display( ASub, "ASub" );
    if( print )
        El::Print( ASub, "ASub" );

    El::MakeUniform( ASub );
    if( display )
        El::Display( ASub, "Scrambled ASub" );
    if( print )
        El::Print( ASub, "Scrambled ASub" );
    El::SetSubmatrix( A, rowInds, colInds, ASub );

    if( display )
        El::Display( A, "Modified Fourier matrix" );
    if( print )
        El::Print( A, "Modified Fourier matrix" );
}

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int n = El::Input("--n","size of matrix",20);
        const El::Int numCorruptRows =
          El::Input("--numCorruptRows","num rows of corrupt submatrix",5);
        const El::Int numCorruptCols =
          El::Input("--numCorruptCols","num cols of corrupt submatrix",5);
        const bool display = El::Input("--display","display matrix?",false);
        const bool print = El::Input("--print","print matrix?",true);
        El::ProcessInput();
        El::PrintInputReport();
        if( numCorruptRows > n || numCorruptCols > n )
            El::LogicError("Corrupt submatrix is too large");

        // Construct a default process grid over COMM_WORLD
        const El::Grid grid( El::mpi::COMM_WORLD );
        
        CorruptFourierSubmatrix<float>
        ( n, numCorruptRows, numCorruptCols, grid, print, display );
        CorruptFourierSubmatrix<double>
        ( n, numCorruptRows, numCorruptCols, grid, print, display );
#ifdef EL_HAVE_QD
        CorruptFourierSubmatrix<El::DoubleDouble>
        ( n, numCorruptRows, numCorruptCols, grid, print, display );
        CorruptFourierSubmatrix<El::QuadDouble>
        ( n, numCorruptRows, numCorruptCols, grid, print, display );
#endif
#ifdef EL_HAVE_QUAD
        CorruptFourierSubmatrix<El::Quad>
        ( n, numCorruptRows, numCorruptCols, grid, print, display );
#endif
#ifdef EL_HAVE_MPC
        CorruptFourierSubmatrix<El::BigFloat>
        ( n, numCorruptRows, numCorruptCols, grid, print, display );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
