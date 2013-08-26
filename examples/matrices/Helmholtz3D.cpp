/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Inverse.hpp"
#include "elemental/matrices/Helmholtz.hpp"
#include "elemental/io.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int nx = Input("--nx","size of x dimension",10);
        const Int ny = Input("--ny","size of y dimension",10);
        const Int nz = Input("--nz","size of z dimension",10);
        const double realShift = Input("--realShift","real part of shift",0.);
        const double imagShift = Input("--imagShift","imag part of shift",0.);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        Complex<double> shift( realShift, imagShift );
        auto H = Helmholtz( DefaultGrid(), nx, ny, nz, shift );
        if( display )
            Display( H, "Helmholtz matrix" );
        if( print )
            Print( H, "Helmholtz matrix:" );

        // (Attempt to) invert the Helmholtz matrix
        Inverse( H );
        if( display )
            Display( H, "Inverse of Helmholtz matrix" );
        if( print )
            Print( H, "Inverse of Helmholtz matrix:" );

        // TODO: Extend to allow for computing SVD of submatrix
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
