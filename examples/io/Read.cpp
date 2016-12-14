/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int m = El::Input("--height","height of matrix",10);
        const El::Int n = El::Input("--width","width of matrix",10);
        const std::string filename =
            El::Input("--filename","filename",std::string(""));
        const bool display = El::Input("--display","display matrix?",true);
        const bool print = El::Input("--print","print matrix?",false);
        El::ProcessInput();
        El::PrintInputReport();

        if( filename == "" )
        {
            El::Output("Please specify a filename to read");
        }
        else
        {
            El::DistMatrix<double> A(m,n);
            El::Read( A, filename );
            if( display )
                El::Display( A, "A (distributed read)" );
            if( print )
                El::Print( A, "A (distributed read)" );
            El::Read( A, filename, El::AUTO, true );
            if( display )
                El::Display( A, "A (sequential read)" );
            if( print )
                El::Print( A, "A (sequential read)" );
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
