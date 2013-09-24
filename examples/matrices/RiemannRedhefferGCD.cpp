/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/matrices/GCDMatrix.hpp"
#include "elemental/matrices/Redheffer.hpp"
#include "elemental/matrices/Riemann.hpp"
#include "elemental/io.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        const Grid& g = DefaultGrid();
        auto A = Riemann<double>( g, n );
        auto B = Redheffer<double>( g, n );
        auto C = GCDMatrix<double>( g, n, n );
        if( display )
        {
            Display( A, "Riemann" );
            Display( B, "Redheffer" );
            Display( C, "GCDMatrix" );
        }
        if( print )
        {
            Print( A, "Riemann matrix:" );
            Print( B, "Redheffer matrix:" );
            Print( C, "GCD matrix:" );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
