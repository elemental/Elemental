/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Real>
void RandomFeasibleLP( El::Int m, El::Int n, El::Int k )
{
    El::Output("Testing with ",El::TypeName<Real>());
    // Create random (primal feasible) inputs for the primal/dual problem
    //    arginf_{x,s} { c^T x | A x = b, G x + s = h, s >= 0 }
    //    argsup_{y,z} { -b^T y - h^T z | A^T y + G^T z + c = 0, z >= 0 }.

    // xFeas and sFeas are only used for problem generation
    El::Matrix<Real> xFeas, sFeas;
    El::Uniform( xFeas, n, 1 ); // Sample over B_1(0)
    El::Uniform( sFeas, k, 1, Real(1), Real(1) ); // Sample over B_1(1)

    El::Matrix<Real> A, G;
    El::Matrix<Real> b, c, h;
    El::Uniform( A, m, n );
    El::Uniform( G, k, n );
    El::Gemv( El::NORMAL, Real(1), A, xFeas, b );
    El::Gemv( El::NORMAL, Real(1), G, xFeas, h );
    h += sFeas;
    El::Uniform( c, n, 1 );

    // Solve the primal/dual Linear Program with the default options
    El::Matrix<Real> x, y, z, s;
    El::Timer timer;
    timer.Start();
    El::LP( A, G, b, c, h, x, y, z, s );
    El::Output("Primal-dual LP took ",timer.Stop()," seconds");

    // Print the primal and dual objective values
    const Real primal = El::Dot(c,x);
    const Real dual = -El::Dot(b,y) - El::Dot(h,z);
    const Real relGap = El::Abs(primal-dual) / El::Max(El::Abs(dual),Real(1));
    El::Output("c^T x = ",primal);
    El::Output("-b^T y - h^T z = ",dual);
    El::Output("|gap| / max( |dual|, 1 ) = ",relGap);

    // Print the relative primal feasibility residual,
    //   || A x - b ||_2 / max( || b ||_2, 1 ).
    El::Matrix<Real> rPrimal;
    El::Gemv( El::NORMAL, Real(1), A, x, rPrimal );
    rPrimal -= b;
    const Real bFrob = El::FrobeniusNorm( b );
    const Real rPrimalFrob = El::FrobeniusNorm( rPrimal );
    const Real primalRelResid = rPrimalFrob / El::Max( bFrob, Real(1) );
    El::Output("|| A x - b ||_2 / || b ||_2 = ",primalRelResid);

    // Print the relative dual feasiability residual,
    //   || A^T y + G^T z + c ||_2 / max( || c ||_2, 1 ).
    El::Matrix<Real> rDual;
    El::Gemv( El::TRANSPOSE, Real(1), A, y, rDual );
    El::Gemv( El::TRANSPOSE, Real(1), G, z, Real(1), rDual );
    rDual += c;
    const Real cFrob = El::FrobeniusNorm( c );
    const Real rDualFrob = El::FrobeniusNorm( rDual );
    const Real dualRelResid = rDualFrob / El::Max( cFrob, Real(1) );
    El::Output
    ("|| A^T y + G^T z + c ||_2 / max( || c ||_2, 1 ) = ",dualRelResid);
    El::Output("");
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int m = El::Input("--m","height of A",70);
        const El::Int n = El::Input("--n","width of A",80);
        const El::Int k = El::Input("--k","height of G",90);
        El::ProcessInput();

        RandomFeasibleLP<float>( m, n, k );
        RandomFeasibleLP<double>( m, n, k );
#ifdef EL_HAVE_QD
        RandomFeasibleLP<El::DoubleDouble>( m, n, k );
        RandomFeasibleLP<El::QuadDouble>( m, n, k );
#endif
#ifdef EL_HAVE_QUAD
        RandomFeasibleLP<El::Quad>( m, n, k );
#endif
#ifdef EL_HAVE_MPC
        RandomFeasibleLP<El::BigFloat>( m, n, k );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
