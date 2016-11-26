/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int m = Input("--m","matrix height",100);
        const Int n = Input("--n","matrix width",100);
        ProcessInput();

        DistMatrix<double> A;
        Uniform( A, m, n );
        Print( A, "A" );
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
