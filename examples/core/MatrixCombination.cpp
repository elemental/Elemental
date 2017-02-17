/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This driver is an example of merging two matrices distributed over
   different sets of processes into a single distributed matrix
   distributed over the union of processes. It is a rewrite of a version
   original due to Yingzhou (Ryan) Li.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename T>
void PerformMerge
( El::Int n,
  const El::Grid& grid,
  const El::Grid& lowerGrid,
  const El::Grid& upperGrid,
  bool print )
{
    if( grid.Rank() == 0 )
        El::Output("Will merge submatrices of type ",El::TypeName<T>());

    El::DistMatrix<T> ALower(lowerGrid), AUpper(upperGrid);
    El::Uniform( ALower, n, n );
    El::Uniform( AUpper, n, n );
    if( print )
    {
        El::Print( ALower, "ALower" );
        El::Print( AUpper, "AUpper" );
    }

    if( grid.Rank() == 0 )
        El::Output("Starting MatrixPartition...");
    El::mpi::Barrier();
    double startTime = El::mpi::Time();

    El::DistMatrix<T> A(grid);
    A.Resize( n, 2*n );
    El::DistMatrix<T> AL(grid), AR(grid);
    El::PartitionRight( A, AL, AR, n );
    AL = ALower;
    AR = AUpper;

    El::mpi::Barrier();
    double stopTime = El::mpi::Time();
    if( grid.Rank() == 0 )
        El::Output("Partition took ",stopTime-startTime," seconds.");

    if( print )
        El::Print( A, "A" );
}

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm;
    const int commSize = El::mpi::Size( comm );
    if( commSize ==1 )
    {
        El::Output("Cannot split one process into two teams");
        return 0;
    }

    try
    {
        const El::Int n = El::Input("--size","size of matrix",100);
        El::Int gridHeight = El::Input("--gridHeight","grid height",0);
        El::Int lowerGridHeight =
          El::Input("--lowerGridHeight","lower grid height",0);
        El::Int upperGridHeight =
          El::Input("--upperGridHeight","upper grid height",0);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = El::Grid::DefaultHeight( commSize );
        El::Grid grid( comm, gridHeight );

        // Split the group for MPI_COMM_WORLD into the lower and upper halves
        // ==================================================================
        El::mpi::Group group;
        El::mpi::CommGroup( comm, group );
        const int lowerGroupSize = commSize/2;
        const int upperGroupSize = commSize - lowerGroupSize;
        std::vector<int> ranks(El::Max(lowerGroupSize,upperGroupSize));
        // Form the lower half
        // -------------------
        El::mpi::Group lowerGroup;
        for( int q=0; q<commSize/2; ++q )
            ranks[q] = q;
        El::mpi::Incl( group, commSize/2, ranks.data(), lowerGroup );
        // Form the upper half
        // -------------------
        El::mpi::Group upperGroup;
        for( int q=0; q<commSize-commSize/2; ++q )
            ranks[q] = q + commSize/2;
        El::mpi::Incl( group, commSize-commSize/2, ranks.data(), upperGroup );

        // Form the lower and upper grids
        // ==============================
        if( lowerGridHeight == 0 )
            lowerGridHeight = El::Grid::DefaultHeight( lowerGroupSize );
        if( upperGridHeight == 0 )
            upperGridHeight = El::Grid::DefaultHeight( upperGroupSize );
        El::Grid lowerGrid( comm, lowerGroup, lowerGridHeight );
        El::Grid upperGrid( comm, upperGroup, upperGridHeight );

        // Perform the merges
        // ==================
        PerformMerge<float>
        ( n, grid, lowerGrid, upperGrid, print );
        PerformMerge<El::Complex<float>>
        ( n, grid, lowerGrid, upperGrid, print );
        PerformMerge<double>
        ( n, grid, lowerGrid, upperGrid, print );
        PerformMerge<El::Complex<double>>
        ( n, grid, lowerGrid, upperGrid, print );
#ifdef EL_HAVE_QD
        PerformMerge<El::DoubleDouble>
        ( n, grid, lowerGrid, upperGrid, print );
        PerformMerge<El::Complex<El::DoubleDouble>>
        ( n, grid, lowerGrid, upperGrid, print );
        PerformMerge<El::QuadDouble>
        ( n, grid, lowerGrid, upperGrid, print );
        PerformMerge<El::Complex<El::QuadDouble>>
        ( n, grid, lowerGrid, upperGrid, print );
#endif
#ifdef EL_HAVE_QUAD
        PerformMerge<El::Quad>
        ( n, grid, lowerGrid, upperGrid, print );
        PerformMerge<El::Complex<El::Quad>>
        ( n, grid, lowerGrid, upperGrid, print );
#endif
#ifdef EL_HAVE_MPC
        PerformMerge<El::BigFloat>
        ( n, grid, lowerGrid, upperGrid, print );
        PerformMerge<El::Complex<El::BigFloat>>
        ( n, grid, lowerGrid, upperGrid, print );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
