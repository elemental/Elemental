/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/matrices/Kahan.hpp"
#include "elemental/io.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of identity matrix",10);
        const double phi = Input("--phi","number in (0,1)",0.2);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        auto A = Kahan( DefaultGrid(), n, phi );
        if( display )
            Display( A, "Kahan" );
        if( print )
            Print( A, "Kahan" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
