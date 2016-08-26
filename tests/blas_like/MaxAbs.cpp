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

// template<typename Real, typename matrix<Real>>
template<typename Real, typename TestMatrix>
void TestMaxAbs(Int m, Int n=1)
{
    DEBUG_ONLY(CallStackEntry cse("TestMaxAbs"))
    TestMatrix A;
    Zeros( A, m, m);
    assert( MaxAbs( A ) == Real(0) );
    Ones( A, m, m);
    assert( MaxAbs( A ) == Real(1) );
    MinIJ( A, m );
    assert( MaxAbs( A ) == Real(m) );
    OneTwoOne( A, m );
    assert( MaxAbs( A ) == Real(2) );
}

void RunTests( Int m)
{
  TestMaxAbs<double,DistMatrix<double>>( m );
  TestMaxAbs<float,DistMatrix<float>>( m );
  TestMaxAbs<double,Matrix<double>>( m );
  TestMaxAbs<float,Matrix<float>>( m );
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    try
    {
	Int m = 1;
	for( Int e = 1; e<4; ++e )
	{
		m *= 10;
		Output(" m = ",m,"\n");
		RunTests( m);
	}
    }
    catch( exception& e ) { ReportException(e); }
    return 0;
}
