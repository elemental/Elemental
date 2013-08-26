/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/matrices/Wilkinson.hpp"
#include "elemental/io.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int k = Input("--order","generate 2k+1 x 2k+1 matrix",5);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        auto W = Wilkinson<double>( DefaultGrid(), k );
        if( display )
        {
            Display( W, "Wilkinson matrix" );
#ifdef HAVE_QT5
            Spy( W, "Wilkinson spy plot" );
#endif
        }
        if( print )
            Print( W, "Wilkinson matrix" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
