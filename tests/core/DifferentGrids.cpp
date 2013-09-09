/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/matrices/Identity.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commSize = mpi::CommSize( comm );
    
    try
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Drop down to a square grid, change the matrix, and redistribute back
        const Int commSqrt = Int(sqrt(double(commSize)));

        std::vector<int> sqrtRanks(commSqrt*commSqrt);
        for( Int i=0; i<commSqrt*commSqrt; ++i )
            sqrtRanks[i] = i;

        mpi::Group group, sqrtGroup;
        
        mpi::CommGroup( comm, group );
        mpi::GroupIncl( group, sqrtRanks.size(), sqrtRanks.data(), sqrtGroup );

        const Grid grid( comm );
        const Grid sqrtGrid( comm, sqrtGroup, grid.Height() );

        DistMatrix<double> A(grid), ASqrt(sqrtGrid);

        Identity( A, m, n );
        if( print )
            Print( A, "A" );

        ASqrt = A;
        if( ASqrt.Participating() )
        {
            if( print )
                Print( ASqrt, "ASqrt := A" );
            Scale( 2., ASqrt );
            if( print )
                Print( ASqrt, "ASqrt := 2 ASqrt" );
        }
        A = ASqrt;
        if( print )
            Print( A, "A := ASqrt" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
