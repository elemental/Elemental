/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/matrices/NormalUniformSpectrum.hpp"
#include "elemental/io.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const double realCenter = Input
            ("--realCenter","real center of uniform eigval distribution",3.);
        const double imagCenter = Input
            ("--imagCenter","imag center of uniform eigval distribution",-4.);
        const double radius = Input
            ("--radius","radius of uniform eigval distribution",2.);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        const Complex<double> center( realCenter, imagCenter );
        auto X = NormalUniformSpectrum( DefaultGrid(), n, center, radius );
        if( display )
            Display( X, "Normal uniform spectrum" );
        if( print )
            Print( X, "X" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
