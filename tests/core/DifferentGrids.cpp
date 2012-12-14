/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <ctime>
#include "elemental.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    
    try
    {
        const int m = Input("--height","height of matrix",100);
        const int n = Input("--width","width of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Drop down to a square grid, change the matrix, and redistribute back
        const int commSqrt = int(sqrt(double(commSize)));

        std::vector<int> sqrtRanks(commSqrt*commSqrt);
        for( int i=0; i<commSqrt*commSqrt; ++i )
            sqrtRanks[i] = i;

        mpi::Group group, sqrtGroup;
        
        mpi::CommGroup( comm, group );
        mpi::GroupIncl( group, sqrtRanks.size(), &sqrtRanks[0], sqrtGroup );

        const Grid grid( comm );
        const Grid sqrtGrid( comm, sqrtGroup );

        DistMatrix<double> A(grid), ASqrt(sqrtGrid);

        Identity( m, n, A );
        if( print )
            A.Print("A");

        ASqrt = A;
        if( print )
            ASqrt.Print("ASqrt := A");

        Scal( 2.0, ASqrt );
        if( print )
            ASqrt.Print("ASqrt := 2 ASqrt");

        A = ASqrt;
        if( print )
            A.Print("A := ASqrt");
    }
    catch( ArgException& e ) { }
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << std::endl;
        std::cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }
    Finalize();
    return 0;
}
