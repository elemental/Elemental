/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/matrices/Cauchy.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int m = Input("--height","height of matrix",10);
        const int n = Input("--width","width of matrix",10);
        const bool print = Input("--print","print matrix?",true);
        ProcessInput();
        PrintInputReport();

        std::vector<double> x( m ), y( n );
        for( int j=0; j<m; ++j )
            x[j] = j;
        for( int j=0; j<n; ++j )
            y[j] = j+m;

        DistMatrix<double> A;
        Cauchy( x, y, A );
        if( print )
            A.Print("Cauchy matrix:");
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
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
