/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_MAKETRIANGULAR_INC
#include ELEM_SETDIAGONAL_INC
#include ELEM_LDL_INC
#include ELEM_WALSH_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

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
            LDLH( W );
            auto d = W.GetDiagonal();
            MakeTriangular( LOWER, W );
            SetDiagonal( W, 1. );

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
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
