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
    const Real epsilon = limits::Epsilon<Real>();

    const Int n = 4;
    Matrix<Real> A;
    Zeros( A, n, n );
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
    auto AOrig( A );

    Matrix<Real> Q;
    Identity( Q, n, n );

    vector<Real> work(n);
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

    // E := I - Q' Q
    Matrix<Real> E;
    Identity( E, n, n );
    Gemm( ADJOINT, NORMAL, Real(-1), Q, Q, Real(1), E );
    const Real orthogErr = OneNorm( E );
    const Real relOrthogErr = orthogErr / epsilon;
    Output("|| I - Q' Q ||_1 / eps = ",relOrthogErr);

    // E := AOrig - Q A Q'
    Matrix<Real> Z;
    Gemm( NORMAL, NORMAL, Real(1), Q, A, Z );
    E = AOrig; 
    Gemm( NORMAL, ADJOINT, Real(-1), Z, Q, Real(1), E );
    const Real oneNormA = OneNorm( AOrig );
    const Real swapErr = OneNorm( E );
    const Real relSwapErr = swapErr / (epsilon*oneNormA);
    Output("|| AOrig - Q A Q' ||_1 / (eps*|| A ||_1) = ",relSwapErr);
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
