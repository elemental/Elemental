/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/matrices/Lehmer.hpp"
#include "elemental/matrices/Parter.hpp"
#include "elemental/matrices/Ris.hpp"
#include "elemental/io.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> A, B, C;
        Lehmer( A, n );
        Parter( B, n );
        Ris( C, n );
        if( display )
        {
            Display( A, "Lehmer" );
            Display( B, "Parter" );
            Display( C, "Ris" );
        }
        if( print )
        {
            A.Print("Lehmer:");
            B.Print("Parter:");
            C.Print("Ris:");
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
