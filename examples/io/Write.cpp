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
        const Int n = Input("--size","size of matrix",100);
        const double omega = Input("--omega","frequency of FoxLi",16*M_PI);
        const std::string basename = 
            Input("--basename","basename of file",std::string(""));
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        if( basename == "" )
        {
            Output("Please specify a basename for writing");
        }
        else
        {
            DistMatrix<Complex<double>> A;
            FoxLi( A, n, omega );
            if( display )
                Display( A, "A" );
            if( print )
                Print( A, "A" );
            Write( A, basename, MATRIX_MARKET );

            DistMatrix<Complex<double>> B;
            Read( B, basename+".mm" );
            if( display )
                Display( B, "B" );
            if( print )
                Print( B, "B" );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
