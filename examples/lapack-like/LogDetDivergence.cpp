/*
   Copyright (c) 2009-2012, Jack Poulson
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
        const int n = Input("--size","size of HPD matrix",100);
        const double lower = Input("--lower","lower bound on spectrum",1.);
        const double upper = Input("--upper","upper bound on spectrum",10.);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> A, B;
        HermitianUniformSpectrum( n, A, lower, upper );
        HermitianUniformSpectrum( n, B, lower, upper );
        if( print )
        {
            A.Print("A");
            B.Print("B");
        }
        const double logDetDiv = LogDetDivergence( LOWER, A, B );
        if( commRank == 0 )
            std::cout << "LogDetDiv(A,B) = " << logDetDiv << std::endl;
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

