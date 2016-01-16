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
#include "El.hpp"
using namespace El;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );
    if( commSize ==1 )
    {
        Output("Cannot split one process into two teams");
        return 0;
    }

    try
    {
        const Int n = Input("--size","size of matrix",100);
        Int gridHeight = Input("--gridHeight","grid height",0);
        Int lowerGridHeight = Input("--lowerGridHeight","lower grid height",0);
        Int upperGridHeight = Input("--upperGridHeight","upper grid height",0);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // If the grid height wasn't specified, then we should attempt to build
        // a nearly-square process grid
        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( commSize );
        Grid grid( comm, gridHeight );

        // Split the group for MPI_COMM_WORLD into the lower and upper halves
        // ==================================================================
        mpi::Group group;
        mpi::CommGroup( comm, group );
        const int lowerGroupSize = commSize/2;
        const int upperGroupSize = commSize - lowerGroupSize;
        std::vector<int> ranks(Max(lowerGroupSize,upperGroupSize));
        // Form the lower half
        // -------------------
        mpi::Group lowerGroup;
        for( int q=0; q<commSize/2; ++q )
            ranks[q] = q;
        mpi::Incl( group, commSize/2, ranks.data(), lowerGroup );
        // Form the upper half
        // -------------------
        mpi::Group upperGroup;
        for( Int q=0; q<commSize-commSize/2; ++q )
            ranks[q] = q + commSize/2;
        mpi::Incl( group, commSize-commSize/2, ranks.data(), upperGroup );
 
        // Form the lower and upper grids
        // ==============================
        if( lowerGridHeight == 0 )
            lowerGridHeight = Grid::FindFactor( lowerGroupSize );
        if( upperGridHeight == 0 )
            upperGridHeight = Grid::FindFactor( upperGroupSize );
        Grid lowerGrid( comm, lowerGroup, lowerGridHeight );
        Grid upperGrid( comm, upperGroup, upperGridHeight );

        // Create the original distributed matrices
        // ========================================
        DistMatrix<double> ALower(lowerGrid), AUpper(upperGrid);
        Uniform( ALower, n, n );
        Uniform( AUpper, n, n );
        if( print )
        {
            Print( ALower, "ALower" );
            Print( AUpper, "AUpper" );
        }

        if( commRank == 0 )
            Output("Starting MatrixPartition...");
        mpi::Barrier();
        double startTime = mpi::Time();

        DistMatrix<double> A(grid);
        A.Resize( n, 2*n );
        DistMatrix<double> AL(grid),AR(grid);
        PartitionRight( A, AL, AR, n );
        AL = ALower;
        AR = AUpper;

        mpi::Barrier();
        double stopTime = mpi::Time();
        if( commRank == 0 )
            Output("Partition took ",stopTime-startTime," seconds.");

        if( print )
            Print( A, "A" );

    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
