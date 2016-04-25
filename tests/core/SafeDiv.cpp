/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        Complex<double> x(Pow(2.,1023.),Pow(2.,-1023.));
        Complex<double> y(Pow(2.,677.),Pow(2.,-677.));
        Complex<double> z = x / y;
        Complex<double> zSafe = SafeDiv( x, y );
        Output("x=",x,", y=",y,", z=",z,", zSafe=",zSafe);
        Output("|x-y*z|/|x|=",Abs(x-y*z)/Abs(x));
        Output("|x-y*zSafe|/|x|=",Abs(x-y*zSafe)/Abs(x));
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
