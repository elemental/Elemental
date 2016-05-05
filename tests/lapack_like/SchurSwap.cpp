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
void TestSwap( Real tau, bool testAccuracy, bool print )
{
    Output("Testing with ",TypeName<Real>());

    Matrix<Real> A;
    Zeros( A, 4, 4 );
    A.Set( 0, 0, Real(7.001) );
    A.Set( 0, 1, Real(-87)   );
    A.Set( 0, 2, Real(39.4)*tau );
    A.Set( 0, 3, Real(22.2)*tau );
    A.Set( 1, 0, Real(5) );
    A.Set( 1, 1, Real(7.001) );
    A.Set( 1, 2, Real(-12.2)*tau );
    A.Set( 1, 3, Real(36.)*tau );
    A.Set( 2, 2, Real(7.01) );
    A.Set( 2, 3, Real(-11.7567) );
    A.Set( 3, 2, Real(37) );
    A.Set( 3, 3, Real(7.01) );
    if( print ) 
        Print( A, "A" );

    Matrix<Real> Q;
    Identity( Q, 4, 4 );

    vector<Real> work(4);
    lapack::AdjacentSchurExchange
    ( A.Height(),
      A.Buffer(), A.LDim(),
      Q.Buffer(), Q.LDim(),
      0, 2, 2,
      work.data(), testAccuracy );
    if( print )
    {
        Print( A, "A" );
        Print( Q, "Q" );
    }

    // TODO: Add accuracy tests; but this has already proven useful
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const double tau = Input("--tau","scaling parameter",1.); 
        const bool testAccuracy = Input("--testAccuracy","test accuracy?",true);
        const bool print = Input("--print","print?",true);
        ProcessInput();

        TestSwap( float(tau), testAccuracy, print );
        TestSwap( double(tau), testAccuracy, print );
#ifdef EL_HAVE_QUAD
        TestSwap( Quad(tau), testAccuracy, print );
#endif
#ifdef EL_HAVE_QD
        TestSwap( DoubleDouble(tau), testAccuracy, print );
        TestSwap( QuadDouble(tau), testAccuracy, print );
#endif
#ifdef EL_HAVE_MPC
        TestSwap( BigFloat(tau), testAccuracy, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
