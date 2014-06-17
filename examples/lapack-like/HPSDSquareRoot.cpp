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
        const Int n = Input("--size","size of HPSD matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> L;
        Uniform( L, n, n );
        MakeTrapezoidal( LOWER, L, -1 );
        DistMatrix<C> A;
        Zeros( A, n, n );
        Herk( LOWER, NORMAL, C(1), L, C(0), A );
        if( print )
            Print( A, "A" );

        // Replace A with its matrix square root
        HPSDSquareRoot( LOWER, A );
        if( print )
        {
            MakeHermitian( LOWER, A );
            Print( A, "sqrt(A)" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
