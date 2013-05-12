/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Norm.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
#include "elemental/graphics.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n = Input("--size","size of Hermitian matrix",10);
        const double lower = Input("--lower","lower bound on spectrum",1.);
        const double upper = Input("--upper","upper bound on spectrum",10.);
        const bool print = Input("--print","print matrix?",true);
#ifdef HAVE_QT5
        const bool display = Input("--display","display matrix?",true);
#endif
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> X;
        HermitianUniformSpectrum( X, n, lower, upper );
        if( print )
            X.Print("X");
#ifdef HAVE_QT5
        if( display )
            Display( X, "Hermitian uniform spectrum" );
#endif
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
