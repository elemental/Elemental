/*
   Copyright (c) 2009-2016, Ryan H. Lewis 
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <random>
#include <stdexcept>

using namespace El;

template<typename T>
void TestMultiply(Int m, Int n=1)
{
    DEBUG_ONLY(CallStackEntry cse("TestMultiply"))
    SparseMatrix<T> A;
    Identity( A, m, m);
    Matrix<T> B( m, n);
    Bernoulli( B, m, n);
    auto G = A.Graph();
    Matrix<T> C(m,n),D(m,n);
    Multiply(NORMAL, T(1), A, B, T(0), C);
    Multiply(NORMAL, T(1), G, B, T(0), D);
    Axpy(T(-1), C, D);
    auto nrm = FrobeniusNorm(D);
    std::cout << "error = " << nrm << std::endl;
    if( nrm > limits::Epsilon<T>()){ RuntimeError("Sparse(I)*x != Graph(I)*x"); }
}

void RunTests( Int m)
{
  TestMultiply<double>( m);
  //List all the types here..
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
