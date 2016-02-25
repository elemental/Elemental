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
        const Int k = Input("--order","generate 2^k x 2^k matrix",4);
        const bool binary = Input("--binary","binary data?",false);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        // Generate a Walsh matrix of order k (a 2^k x 2^k matrix)
        DistMatrix<double> W;
        Walsh( W, k, binary );
        if( display )
            Display( W, "Walsh matrix" );
        if( print )
            Print( W, "W(2^k)");

        if( !binary )
        {
            LDL( W, true );
            auto d = GetDiagonal(W);
            MakeTrapezoidal( LOWER, W );
            FillDiagonal( W, 1. );

            if( display )
            {
                Display( W, "Lower factor" );
                Display( d, "Diagonal factor" );
            }
            if( print )
            {
                Print( W, "L" ); 
                Print( d, "d" ); 
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
