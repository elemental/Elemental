/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename Real>
void Test()
{
    Output("Testing with ",TypeName<Complex<Real>>());
    PushIndent();
    Complex<Real> x(Pow(Real(2),Real(1023)),Pow(Real(2),Real(-1023)));
    Complex<Real> y(Pow(Real(2),Real(677)),Pow(Real(2),Real(-677)));
    Complex<Real> z = x / y;
    Complex<Real> zSafe = SafeDiv( x, y );
    Output("x=",x,", y=",y,", z=",z,", zSafe=",zSafe);
    Output("|x-y*z|/|x|=",Abs(x-y*z)/Abs(x));
    Output("|x-y*zSafe|/|x|=",Abs(x-y*zSafe)/Abs(x));
    PopIndent();
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        Test<double>();
#ifdef EL_HAVE_QUAD
        Test<Quad>();
#endif
#ifdef EL_HAVE_MPC
        Test<BigFloat>();
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
