/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","matrix height",100);
        const Int n = Input("--width","matrix width",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Complex<double>> A, Q, P;
        Uniform( A, m, n );

        // Compute the polar decomp of A (but do not overwrite A)
        Q = A;
        Polar( Q, P );

        if( print )
        {
            Print( A, "A" );
            Print( Q, "Q" );
            Print( P, "P" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
