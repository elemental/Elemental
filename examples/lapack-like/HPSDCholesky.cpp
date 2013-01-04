/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try 
    {
        const int n = Input("--size","size of HPSD matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        Grid g( comm );
        DistMatrix<C> L(g), A(g);
        Uniform( n, n, L );
        MakeTrapezoidal( LEFT, LOWER, -1, L );
        Zeros( n, n, A );
        Herk( LOWER, NORMAL, C(1), L, C(0), A );

        if( print )
            A.Print("A");

        // Replace A with its Cholesky factor
        HPSDCholesky( LOWER, A );

        if( print )
        {
            MakeTriangular( LOWER, A );
            A.Print("chol(A)");
        }
    }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught exception with message: "
           << e.what() << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}

