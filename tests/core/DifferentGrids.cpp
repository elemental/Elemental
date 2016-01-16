/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

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
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commSize = mpi::Size( comm );
    
    try
    {
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const bool colMajorSqrt = Input("--colMajorSqrt","colMajor sqrt?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const GridOrder orderSqrt = ( colMajorSqrt ? COLUMN_MAJOR : ROW_MAJOR );

        // Drop down to a square grid, change the matrix, and redistribute back
        const Int commSqrt = Int(sqrt(double(commSize)));

        std::vector<int> sqrtRanks(commSqrt*commSqrt);
        for( Int i=0; i<commSqrt*commSqrt; ++i )
            sqrtRanks[i] = i;

        mpi::Group group, sqrtGroup;
        
        mpi::CommGroup( comm, group );
        mpi::Incl( group, sqrtRanks.size(), sqrtRanks.data(), sqrtGroup );

        const Grid grid( comm, order );
        const Grid sqrtGrid( comm, sqrtGroup, commSqrt, orderSqrt );

        DistMatrix<double> A(grid), ASqrt(sqrtGrid);

        Identity( A, m, n );
        if( print )
            Print( A, "A" );

        ASqrt = A;
        if( ASqrt.Participating() )
        {
            if( print )
                Print( ASqrt, "ASqrt := A" );
            ASqrt *= 2;
            if( print )
                Print( ASqrt, "ASqrt := 2 ASqrt" );
        }
        A = ASqrt;
        if( print )
            Print( A, "A := ASqrt" );

        const Grid newGrid( comm, order );
        A.SetGrid( newGrid );
        if( print )
            Print( A, "A after changing grid" );
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
