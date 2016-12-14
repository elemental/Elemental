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
        const El::Int n = El::Input("--size","size of matrix",100);
        const double omega = El::Input("--omega","frequency of FoxLi",16*M_PI);
        const std::string basename =
            El::Input("--basename","basename of file",std::string(""));
        const bool display = El::Input("--display","display matrix?",true);
        const bool print = El::Input("--print","print matrix?",false);
        El::ProcessInput();
        El::PrintInputReport();

        if( basename == "" )
        {
            El::Output("Please specify a basename for writing");
        }
        else
        {
            El::DistMatrix<El::Complex<double>> A;
            El::FoxLi( A, n, omega );
            if( display )
                El::Display( A, "A" );
            if( print )
                El::Print( A, "A" );
            El::Write( A, basename, El::MATRIX_MARKET );

            El::DistMatrix<El::Complex<double>> B;
            El::Read( B, basename+".mm" );
            if( display )
                El::Display( B, "B" );
            if( print )
                El::Print( B, "B" );
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
