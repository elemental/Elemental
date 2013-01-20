/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Norm.hpp"
#include "elemental/lapack-like/Pseudoinverse.hpp"
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
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        Grid g( comm );
        DistMatrix<C> A( g );
        Uniform( m, n, A );

        // Compute the pseudoinverseof A (but do not overwrite A)
        DistMatrix<C> pinvA( A );
        Pseudoinverse( pinvA );

        if( print )
        {
            A.Print("A");
            pinvA.Print("pinv(A)");
        }

        const R frobOfA = Norm( A, FROBENIUS_NORM );
        const R frobOfPinvA = Norm( pinvA, FROBENIUS_NORM );

        if( commRank == 0 )
        {
            cout << "||   A   ||_F =  " << frobOfA << "\n"
                 << "||pinv(A)||_F =  " << frobOfPinvA << "\n"
                 << endl;
        }
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
