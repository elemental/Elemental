/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/lapack-like/LDL.hpp"
#include "elemental/matrices/Walsh.hpp"
#include "elemental/io.hpp"
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
        auto W = Walsh<double>( DefaultGrid(), k, binary );
        if( display )
            Display( W, "Walsh matrix" );
        if( print )
            Print( W, "W(2^k)");

        if( !binary )
        {
            DistMatrix<double,MC,STAR> d;
            LDLH( W, d );
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
