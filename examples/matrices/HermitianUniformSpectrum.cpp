/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Norm.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
#include "elemental/io.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of Hermitian matrix",10);
        const double lower = Input("--lower","lower bound on spectrum",1.);
        const double upper = Input("--upper","upper bound on spectrum",10.);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        const Grid& g = DefaultGrid();
        auto X = HermitianUniformSpectrum<double>( g, n, lower, upper );
        if( display )
            Display( X, "Hermitian uniform spectrum" );
        if( print )
            Print( X, "X" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
