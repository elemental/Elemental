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
        const int n = Input("--size","size of matrix",10);
        const double realCenter = Input
            ("--realCenter","real center of uniform eigval distribution",3.);
        const double imagCenter = Input
            ("--imagCenter","imag center of uniform eigval distribution",-4.);
        const double radius = Input
            ("--radius","radius of uniform eigval distribution",2.);
        const bool print = Input("--print","print matrix?",true);
        ProcessInput();
        PrintInputReport();

        const Complex<double> center( realCenter, imagCenter );
        DistMatrix<Complex<double> > X;
        NormalUniformSpectrum( n, X, center, radius );
        if( print )
            X.Print("X");
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
