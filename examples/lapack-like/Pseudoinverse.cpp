/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"

using namespace std;
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );

        // Compute the pseudoinverseof A (but do not overwrite A)
        DistMatrix<C> pinvA( A );
        Pseudoinverse( pinvA );
        if( print )
        {
            Print( A, "A" );
            Print( pinvA, "pinv(A)" );
        }

        const Real frobOfA = FrobeniusNorm( A );
        const Real frobOfPinvA = FrobeniusNorm( pinvA );

        if( mpi::WorldRank() == 0 )
        {
            cout << "||   A   ||_F =  " << frobOfA << "\n"
                 << "||pinv(A)||_F =  " << frobOfPinvA << "\n"
                 << endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
