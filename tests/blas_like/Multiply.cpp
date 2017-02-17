/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2016, Ryan H. Lewis
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
    EL_DEBUG_CSE
    typedef Base<T> Real;
    Output("Testing with ",TypeName<T>());

    SparseMatrix<T> A;
    Identity( A, m, m );
    auto G = A.Graph();

    Matrix<T> B;
    Bernoulli( B, m, n );

    Matrix<T> C, D;
    Zeros( C, m, n );
    Zeros( D, m, n );
    Multiply( NORMAL, T(1), A, B, T(0), C );
    Multiply( NORMAL, T(1), G, B, T(0), D );
    Axpy( T(-1), C, D );
    auto DFrob = FrobeniusNorm(D);
    if( DFrob > limits::Epsilon<Real>() )
    {
        Output("|| A B - G B ||_F = ",DFrob);
        RuntimeError("Sparse(I)*x != Graph(I)*x");
    }
    else
        Output("Test passed");
}

void RunTests( Int m )
{
    PushIndent();
    TestMultiply<float>(m);
    TestMultiply<Complex<float>>(m);
    TestMultiply<double>(m);
    TestMultiply<Complex<double>>(m);
#ifdef EL_HAVE_QD
    TestMultiply<DoubleDouble>(m);
    TestMultiply<Complex<DoubleDouble>>(m);
    TestMultiply<QuadDouble>(m);
    TestMultiply<Complex<QuadDouble>>(m);
#endif
#ifdef EL_HAVE_QUAD
    TestMultiply<Quad>(m);
    TestMultiply<Complex<Quad>>(m);
#endif
#ifdef EL_HAVE_MPC
    TestMultiply<BigFloat>(m);
    TestMultiply<Complex<BigFloat>>(m);
#endif
    PopIndent();
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    try
    {
        Int m = 1;
        for( Int e=1; e<4; ++e )
        {    
            m *= 10;    
            Output("Testing with matrix height of ",m);
            RunTests(m);
        }
    }
    catch( exception& e ) { ReportException(e); }
    return 0;
}
