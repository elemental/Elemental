/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_INVERSE_INC
#include ELEM_HELMHOLTZ_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--n","size of matrix",100);
        const double realShift = Input("--realShift","real part of shift",0.);
        const double imagShift = Input("--imagShift","imag part of shift",0.);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        const bool write = Input("--write","write matrix?",false);
        const Int formatInt = Input("--format","write format",2);
        const Int colorMapInt = Input("--colorMap","color map",0);
        ProcessInput();
        PrintInputReport();

        if( formatInt < 1 || formatInt >= FileFormat_MAX )
            LogicError("Invalid file format integer, should be in [1,",
                       FileFormat_MAX,")");

        FileFormat format = static_cast<FileFormat>(formatInt);
        ColorMap colorMap = static_cast<ColorMap>(colorMapInt);
        SetColorMap( colorMap );

        Complex<double> shift( realShift, imagShift );
        DistMatrix<Complex<double>> H;
        Helmholtz( H, n, shift );
        if( display )
            Display( H, "Helmholtz matrix" );
        if( print )
            Print( H, "Helmholtz matrix:" );
        if( write )
            Write( H, "H", format );

        // (Attempt to) invert the Helmholtz matrix
        Inverse( H );
        if( display )
            Display( H, "Inverse of Helmholtz matrix" );
        if( print )
            Print( H, "Inverse of Helmholtz matrix:" );
        if( write )
            Write( H, "invH", format );

        // TODO: Extend to allow for computing SVD of submatrix
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
