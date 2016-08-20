/*
   Copyright (c) 2009-2016, Jack Poulson
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
void TestMultiply(Int m, Int n=100)
{
    DEBUG_ONLY(CallStackEntry cse("TestMultiply"))
    SparseMatrix<T> A;
    Identity( A, m, m);
    Matrix<T> B;
    Bernoulli( B, m, n);
    auto G = A.Graph();
    Matrix<T> C,D;
    Reshape(m, n, C);
    Reshape(m, n, D);
    Multiply(NORMAL, 1.0, A, B, 0.0, C);
    Multiply(NORMAL, 1.0, G, B, 0.0, D);
    Axpy(-1.0, C, D);
    auto nrm = FrobeniusNorm(D);
    if( nrm > 1e-16){ throw std::runtime_error("Sparse(I)*x != Graph(I)*x"); }
}

template<typename T>
void TestMultiplyDiagonal(Int m, Int n=100)
{
    DEBUG_ONLY(CallStackEntry cse("TestMultiplyDiagonal"))
    std::random_device rd;
    std::mt19937 gen(rd());
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<> d(5,2);
    SparseMatrix<T> A;
    Laplacian( A, m);
    Matrix<T> B;
    Bernoulli( B, m, n);
    auto G = A.Graph();
    Matrix<T> C,D;
    Reshape(m, n, C);
    Reshape(m, n, D);
    Multiply(NORMAL, 1.0, A, B, 0.0, C);
    Multiply(NORMAL, 1.0, G, B, 0.0, D);
    Axpy(-1.0, C, D);
    auto nrm = FrobeniusNorm(D);
    if( nrm > 1e-16){ throw std::runtime_error("Sparse(J)*x != Graph(J)*x"); }
}

void RunTests( Int m){
  TestMultiply<float>( m);
  TestMultiply<double>( m);
  TestMultiplyDiagonal<float>( m);
  TestMultiplyDiagonal<double>( m);	
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
	Int m = 1;
	for(int e = 1; e < 4; ++e){	
		m *= 10;
		RunTests( m);
	}
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
