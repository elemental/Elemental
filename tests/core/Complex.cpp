/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try 
    {
        Complex<double> w( 0, 0 );
        if( commRank == 0 )
            std::cout << "sqrt of " << w << " is " << Sqrt(w) << std::endl;

        double maxLocalError = 0;
        double maxLocalRelError = 0;
        const int numTests = 1000;
        for( int j=0; j<numTests; ++j )
        {
            w = SampleUnitBall<Complex<double> >();
            Complex<double> sqrtW = Sqrt(w);
            const double error = Abs(sqrtW*sqrtW-w);
            const double relError = error/Abs(w);

            maxLocalError = std::max( maxLocalError, error );
            maxLocalRelError = std::max( maxLocalRelError, relError );
        }

        double maxError, maxRelError;
        mpi::Reduce( &maxLocalError, &maxError, 1, mpi::SUM, 0, comm );
        mpi::Reduce( &maxLocalRelError, &maxRelError, 1, mpi::SUM, 0, comm );
        if( commRank == 0 )
            std::cout << "Maximum error and relative error from " << numTests 
                      << " tests was " << maxError << " and " << maxRelError
                      << std::endl;
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

