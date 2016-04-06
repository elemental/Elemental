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

#ifdef EL_HAVE_QD
    try
    {
        DoubleDouble a( 17.1374 );
        Output("a=",a,", round(a)=",static_cast<long long int>(a));

        a = 1e18;
        Output("a=",a,", round(a)=",static_cast<long long int>(a));

        // Test an overflow result
        a = 1e23;
        Output("a=",a,", round(a)=",static_cast<long long int>(a));
        

        QuadDouble b( 17.1374 );
        Output("b=",b,", round(b)=",static_cast<long long int>(b));

        b = 1e18;
        Output("b=",b,", round(b)=",static_cast<long long int>(b));

        // Test an overflow result
        b = 1e23;
        Output("b=",b,", round(b)=",static_cast<long long int>(b));
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
