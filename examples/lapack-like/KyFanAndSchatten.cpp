/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Norm/Entrywise.hpp"
#include "elemental/lapack-like/Norm/KyFan.hpp"
#include "elemental/lapack-like/Norm/Schatten.hpp"
#include "elemental/matrices/Uniform.hpp"
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
        const int m = Input("--height","height of matrix",100);
        const int n = Input("--width","width of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const int k = Input("--k","index of KyFan norm",10);
        const double p = Input("--p","power of Schatten norm",2);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );

        Grid g( comm );
        DistMatrix<C> A( g );
        Uniform( m, n, A );

        if( print )
            A.Print("A");

        const double kyFanNorm = KyFanNorm( A, k );
        const double schattenNorm = SchattenNorm( A, p );
        const double entrywiseNorm = EntrywiseNorm( A, p );
        if( commRank == 0 )
            cout << "|| A ||_K(p)   = " << kyFanNorm << "\n"
                 << "|| A ||_S(p)   = " << schattenNorm << "\n"
                 << "|| vec(A) ||_p = " << entrywiseNorm << std::endl;
    }
    catch( ArgException& e )
    {
        // There is nothing to do
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
