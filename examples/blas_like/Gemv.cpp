/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Field>
void TestGemv
( El::Int height,
  El::Int width,
  El::Orientation orientation,
  const El::Grid& grid,
  bool print )
{
    El::DistMatrix<Field> A(grid);
    El::Uniform( A, height, width );

    // Draw the entries of the original x and y from uniform distributions
    // over the complex unit ball
    El::DistMatrix<Field,El::VC,El::STAR> x(grid), y(grid);
    if( orientation == El::NORMAL )
    {
        El::Uniform( x, width, 1 );
        El::Uniform( y, height, 1 );
    }
    else
    {
        El::Uniform( x, height, 1 );
        El::Uniform( y, width, 1 );
    }

    if( print )
    {
        El::Print( A, "A" );
        El::Print( x, "x" );
        El::Print( y, "y" );
    }

    // Run the matrix-vector product
    if( grid.Rank() == 0 )
        El::Output("Starting Gemv (with Field=",El::TypeName<Field>(),")");
    El::Timer gemvElem;
    gemvElem.Start();
    // Form y := 3 A x + 4 y
    El::Gemv( orientation, Field(3), A, x, Field(4), y );
    gemvElem.Stop();
    if( grid.Rank() == 0 )
        El::Output("  Time: ",gemvElem.Total());

    if( print )
    {
        if( orientation == El::NORMAL )
            El::Print( y, "y := 3 A x + 4 y" );
        else
            El::Print( y, "y := 3 A^H x + 4 y" );
    }
}

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int height = El::Input("--height","height of matrix",100);
        const El::Int width = El::Input("--width","width of matrix",100);
        const El::Int blocksize =
          El::Input("--blocksize","algorithmic blocksize",96);
        const bool adjoint = El::Input("--adjoint","apply adjoint?",false);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::SetBlocksize( blocksize );
        const El::Orientation orientation = adjoint ? El::ADJOINT : El::NORMAL;

        // Construct a default grid over the entire set of processes
        const El::Grid grid( El::mpi::COMM_WORLD );

        TestGemv<float>
        ( height, width, orientation, grid, print );
        TestGemv<El::Complex<float>>
        ( height, width, orientation, grid, print );
        TestGemv<double>
        ( height, width, orientation, grid, print );
        TestGemv<El::Complex<double>>
        ( height, width, orientation, grid, print );
#ifdef EL_HAVE_QD
        TestGemv<El::DoubleDouble>
        ( height, width, orientation, grid, print );
        TestGemv<El::Complex<El::DoubleDouble>>
        ( height, width, orientation, grid, print );

        TestGemv<El::QuadDouble>
        ( height, width, orientation, grid, print );
        TestGemv<El::Complex<El::QuadDouble>>
        ( height, width, orientation, grid, print );
#endif
#ifdef EL_HAVE_QUAD
        TestGemv<El::Quad>
        ( height, width, orientation, grid, print );
        TestGemv<El::Complex<El::Quad>>
        ( height, width, orientation, grid, print );
#endif
#ifdef EL_HAVE_MPC
        TestGemv<El::BigFloat>
        ( height, width, orientation, grid, print );
        TestGemv<El::Complex<El::BigFloat>>
        ( height, width, orientation, grid, print );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
