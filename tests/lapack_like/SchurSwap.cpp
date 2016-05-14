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
void TestSwap( const Matrix<Real>& AOrig, bool testAccuracy, bool print )
{
    const Real epsilon = limits::Epsilon<Real>();

    auto A( AOrig );
    const Int n = A.Height();
    if( n != 4 )
        LogicError("Currently assuming n=4");
    if( print ) 
        Print( A, "A" );

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

template<typename Real>
void Test( Real tau, bool testAccuracy, bool print )
{
    Output("Testing with ",TypeName<Real>());
    PushIndent();
    const Int n = 4;

    // Try several examples from Bai and Demmel's paper
    Output("Trying 'tau' example with tau=",tau);
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
    TestSwap( A, testAccuracy, print );

    // Test 1 from Table 1
    Output("Trying Test 1 from Table 1");
    Zeros( A, n, n );    
    A.Set( 0, 0, Real(2) );
    A.Set( 0, 1, Real(-87) );
    A.Set( 0, 2, Real(-20000) );
    A.Set( 0, 3, Real(10000) );
    A.Set( 1, 0, Real(5) );
    A.Set( 1, 1, Real(2) );
    A.Set( 1, 2, Real(-20000) );
    A.Set( 1, 3, Real(-10000) );
    A.Set( 2, 2, Real(1) );
    A.Set( 2, 3, Real(-11) );
    A.Set( 3, 2, Real(37) );
    A.Set( 3, 3, Real(1) );
    TestSwap( A, testAccuracy, print );

    // Test 2 from Table 1
    Output("Trying Test 2 from Table 1");
    Zeros( A, n, n );
    A.Set( 0, 0, Real(1) );
    A.Set( 0, 1, Real(-3) );
    A.Set( 0, 2, Real(3576) );
    A.Set( 0, 3, Real(4888) );
    A.Set( 1, 0, Real(1) );
    A.Set( 1, 1, Real(1) );
    A.Set( 1, 2, Real(-88) );
    A.Set( 1, 3, Real(-1440) );
    A.Set( 2, 2, Real(1.001) );
    A.Set( 2, 3, Real(-3) );
    A.Set( 3, 2, Real(1.001) );
    A.Set( 3, 3, Real(1.001) );
    TestSwap( A, testAccuracy, print );

    // Test 3 from Table 1
    Output("Trying Test 3 from Table 1");
    Zeros( A, n, n );
    A.Set( 0, 0, Real(1) );
    A.Set( 0, 1, Real(-100) );
    A.Set( 0, 2, Real(400) );
    A.Set( 0, 3, Real(-1000) );
    A.Set( 1, 0, Real(0.01) );
    A.Set( 1, 1, Real(1) );
    A.Set( 1, 2, Real(1200) );
    A.Set( 1, 3, Real(-10) );
    A.Set( 2, 2, Real(1.001) );
    A.Set( 2, 3, Real(-0.01) );
    A.Set( 3, 2, Real(100) );
    A.Set( 3, 3, Real(1.001) );
    TestSwap( A, testAccuracy, print );

    // Test 4 from Table 1
    Output("Trying Test 4 from Table 1");
    Zeros( A, n, n );
    A.Set( 0, 0, Real(1) );
    A.Set( 0, 1, Real(-3) );
    A.Set( 0, 2, Real(3) );
    A.Set( 0, 3, Real(2) );
    A.Set( 1, 0, Real(1) );
    A.Set( 1, 1, Real(9) );
    A.Set( 1, 2, Real(0) );
    A.Set( 2, 2, Real(1) );
    A.Set( 2, 3, Real(-3) );
    A.Set( 3, 2, Real(1) );
    A.Set( 3, 3, Real(1) );
    TestSwap( A, testAccuracy, print );
    
    PopIndent();
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const double tau = Input("--tau","scaling parameter",1.); 
        const bool testAccuracy = Input("--testAccuracy","test accuracy?",true);
        const bool print = Input("--print","print?",false);
        ProcessInput();

        Test( float(tau), testAccuracy, print );
        Test( double(tau), testAccuracy, print );
#ifdef EL_HAVE_QUAD
        Test( Quad(tau), testAccuracy, print );
#endif
#ifdef EL_HAVE_QD
        Test( DoubleDouble(tau), testAccuracy, print );
        Test( QuadDouble(tau), testAccuracy, print );
#endif
#ifdef EL_HAVE_MPC
        Test( BigFloat(tau), testAccuracy, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
