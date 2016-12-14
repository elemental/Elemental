/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <cassert>
using namespace El;

template<typename Real, typename TestMatrix>
void TestMaxAbs( Int m )
{
    EL_DEBUG_ONLY(CallStackEntry cse("TestMaxAbs"))
    TestMatrix A;
    Zeros( A, m, m);
    assert( MaxAbs( A ) == Real(0) );
    Ones( A, m, m);
    assert( MaxAbs( A ) == Real(1) );
    Scale( Real(5), A );
    assert( MaxAbs( A ) == Real(5) );
    Scale( Real(-5), A );
    assert( MaxAbs( A ) == Real(25) );
    MinIJ( A, m );
    assert( MaxAbs( A ) == Real(m) );
    OneTwoOne( A, m );
    assert( MaxAbs( A ) == Real(2) );
}

template<typename F, typename TestMatrix>
void TestComplexMaxAbs( Int m )
{
    EL_DEBUG_ONLY(CallStackEntry cse("TestComplexMaxAbs"))
    TestMatrix A;
    Ones( A, m, m );
    assert( ( MaxAbs( A ) - Base<F>( 2.828427125 ) ) < Base<F>(1e-8) );
    Scale( Base<F>( -1 ), A );
    assert( ( MaxAbs( A ) - Base<F>( 2.828427125 ) ) < Base<F>(1e-8) );
}

template<typename Real, typename TestMatrix>
void TestSymmetricMaxAbs( Int m )
{
    EL_DEBUG_ONLY(CallStackEntry cse("TestSymmetricMaxAbs"))
    TestMatrix A;
    Zeros( A, m, m);
    assert( SymmetricMaxAbs( UpperOrLower::LOWER, A ) == Real(0) );
    Ones( A, m, m);
    assert( SymmetricMaxAbs( UpperOrLower::LOWER, A ) == Real(1) );
    Scale( Real(5), A );
    assert( SymmetricMaxAbs( UpperOrLower::LOWER, A ) == Real(5) );
    Scale( Real(-5), A );
    assert( SymmetricMaxAbs( UpperOrLower::LOWER, A ) == Real(25) );
    MinIJ( A, m );
    assert( SymmetricMaxAbs( UpperOrLower::LOWER, A ) == Real(m) );
    OneTwoOne( A, m );
    assert( SymmetricMaxAbs( UpperOrLower::LOWER, A ) == Real(2) );
    Zeros( A, m, m);
    assert( SymmetricMaxAbs( UpperOrLower::UPPER, A ) == Real(0) );
    Ones( A, m, m);
    assert( SymmetricMaxAbs( UpperOrLower::UPPER, A ) == Real(1) );
    Scale( Real(5), A );
    assert( SymmetricMaxAbs( UpperOrLower::UPPER, A ) == Real(5) );
    Scale( Real(-5), A );
    assert( SymmetricMaxAbs( UpperOrLower::UPPER, A ) == Real(25) );
    MinIJ( A, m );
    assert( SymmetricMaxAbs( UpperOrLower::UPPER, A ) == Real(m) );
    OneTwoOne( A, m );
    assert( SymmetricMaxAbs( UpperOrLower::UPPER, A ) == Real(2) );
}

void RunTests( Int m)
{
    TestMaxAbs<double,DistMatrix<double>>( m );
    TestMaxAbs<float,DistMatrix<float>>( m );
    TestMaxAbs<double,Matrix<double>>( m );
    TestMaxAbs<float,Matrix<float>>( m );

    TestComplexMaxAbs<Complex<double>, DistMatrix<Complex<double>>>( m );
    TestComplexMaxAbs<Complex<double>, Matrix<Complex<double>>>( m );
    TestComplexMaxAbs<Complex<float>, DistMatrix<Complex<float>>>( m );
    TestComplexMaxAbs<Complex<float>, Matrix<Complex<float>>>( m );

    TestSymmetricMaxAbs<double,DistMatrix<double>>( m );
    TestSymmetricMaxAbs<float,DistMatrix<float>>( m );
    TestSymmetricMaxAbs<double,Matrix<double>>( m );
    TestSymmetricMaxAbs<float,Matrix<float>>( m );
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    try
    {
        Int m = 1;
        for( Int e = 1; e<4; ++e )
        {
            m *= 10;
            Output(" m = ",m,"\n");
            RunTests( m );
        }
    }
    catch( exception& e ) { ReportException(e); }
    return 0;
}
