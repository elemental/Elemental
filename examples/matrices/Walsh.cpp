/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/lapack-like/LDL.hpp"
#include "elemental/matrices/Walsh.hpp"
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
        const int k = Input("--order","generate 2^k x 2^k matrix",4);
        const bool binary = Input("--binary","binary data?",false);
        const bool print = Input("--print","print matrix?",true);
#ifdef HAVE_QT5
        const bool display = Input("--display","display matrix?",true);
#endif
        ProcessInput();
        PrintInputReport();

        // Generate a Walsh matrix of order k (a 2^k x 2^k matrix)
        DistMatrix<double> W;
        Walsh( W, k, binary );
        if( print )
            W.Print("W(2^k)");
#ifdef HAVE_QT5
        if( display )
            Display( W, "Walsh matrix" );
#endif

        if( !binary )
        {
            DistMatrix<double,MC,STAR> d;
            LDLH( W, d );
            MakeTriangular( LOWER, W );
            SetDiagonal( LEFT, 0, W, 1. );

            if( print )
            {
                W.Print("L"); 
                d.Print("d"); 
            }
#ifdef HAVE_QT5
            if( display )
            {
                Display( W, "Lower factor" );
                Display( d, "Diagonal factor" );
            }
#endif
        }
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
